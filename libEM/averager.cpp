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
const string VarianceAverager::NAME = "variance";
const string TomoAverager::NAME = "mean.tomo";
const string MinMaxAverager::NAME = "minmax";
const string AbsMaxMinAverager::NAME = "absmaxmin";
const string IterationAverager::NAME = "iteration";
const string CtfCWautoAverager::NAME = "ctfw.auto";
const string CtfCAutoAverager::NAME = "ctf.auto";
const string CtfWtAverager::NAME = "ctf.weight";
const string CtfWtFiltAverager::NAME = "ctf.weight.autofilt";
const string FourierWeightAverager::NAME = "weightedfourier";
const string LocalWeightAverager::NAME = "localweight";

template <> Factory < Averager >::Factory()
{
	force_add<ImageAverager>();
	force_add<VarianceAverager>();
	force_add<MinMaxAverager>();
	force_add<AbsMaxMinAverager>();
	force_add<LocalWeightAverager>();
	force_add<IterationAverager>();
	force_add<CtfCWautoAverager>();
	force_add<CtfCAutoAverager>();
	force_add<CtfWtAverager>();
	force_add<CtfWtFiltAverager>();
	force_add<TomoAverager>();
	force_add<FourierWeightAverager>();
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

// pure virtual causing boost issues here, so this is just a dummy method
void Averager::add_image(EMData* image)
{
return;
}

void Averager::add_image_list(const vector<EMData*> & image_list)
{
	for (size_t i = 0; i < image_list.size(); i++) {
		add_image(image_list[i]);
	}
}

TomoAverager::TomoAverager()
	: norm_image(0),nimg(0),overlap(0)
{

}

void TomoAverager::add_image(EMData * image)
{
	if (!image) {
		return;
	}

	if (!image->is_complex()) {
		image=image->do_fft();
		image->set_attr("free_me",(int)1);
	}
		
	if (result!=0 && !EMUtil::is_same_size(image, result)) {
		LOGERR("%s Averager can only process same-size Images",
			   get_name().c_str());
		return;
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
//	size_t image_size = (size_t)nx * ny * nz;

	if (norm_image == 0) {
//		printf("init average %d %d %d",nx,ny,nz);
		result = image->copy_head();
		result->to_zero();
		
		norm_image = image->copy_head();
		norm_image->to_zero();
		
		thresh_sigma = (float)params.set_default("thresh_sigma", 0.5);
		overlap=0.0f;
		nimg=0;
	}

	float *result_data = result->get_data();
	float *norm_data = norm_image->get_data();
	float *data = image->get_data();
	
	vector<float> threshv;
	threshv=image->calc_radial_dist(nx/2,0,1,4);
	for (int i=0; i<nx/2; i++) threshv[i]*=threshv[i]*thresh_sigma;  // The value here is amplitude, we square to make comparison less expensive
	
	size_t j=0;
	// Add any values above threshold to the result image, and add 1 to the corresponding pixels in the norm image
	int k=0;
	for (int z=0; z<nz; z++) {
		for (int y=0; y<ny; y++) {
			for (int x=0; x<nx; x+=2, j+=2) {
				float rf=Util::hypot3(x/2,y<ny/2?y:ny-y,z<nz/2?z:nz-z);	// origin at 0,0; periodic
				int r=int(rf);
				if (r>ny/2) continue;
				float f=data[j];	// real
				float g=data[j+1];	// imag
				float inten=f*f+g*g;
				
				if (inten<threshv[r]) continue;
				
				k+=1;
				result_data[j]  +=f;
				result_data[j+1]+=g;
				
				norm_data[j]  +=1.0;
				norm_data[j+1]+=1.0;
			}
		}
	}
//	printf("%d %d\n",k,nx*ny*nz);
	overlap+=(float)k/(nx*ny*nz);
	nimg++;
	
	if (image->has_attr("free_me")) delete image;
}

EMData * TomoAverager::finish()
{
	if (norm_image==0 || result==0 || nimg==0) return NULL;
	
	int nx = result->get_xsize();
	int ny = result->get_ysize();
	int nz = result->get_zsize();
	size_t image_size = (size_t)nx * ny * nz;
	
	float *result_data = result->get_data();
	float *norm_data = norm_image->get_data();
	
//	printf("finish average %d %d %d",nx,ny,nz);
	// normalize the average
	for (size_t j = 0; j < image_size; j++) {
		if (norm_data[j]==0.0) result_data[j]=0.0;
		else result_data[j]/=norm_data[j];
	}
	
	norm_image->update();
	result->update();
	
	EMData *ret = result->do_ift();
	ret->set_attr("ptcl_repr",norm_image->get_attr("maximum"));
	ret->set_attr("mean_coverage",(float)(overlap/nimg));
	if ((int)params.set_default("save_norm", 0)) 
		norm_image->write_image("norm.hdf");
	
	delete result;
	delete norm_image;
	result=0;
	norm_image=0;
	
	return ret;
}


ImageAverager::ImageAverager()
	: sigma_image(0), ignore0(0), normimage(0), freenorm(0), nimg(0)
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
	size_t image_size = (size_t)nx * ny * nz;

	if (nimg == 1) {
		result = image->copy_head();
		result->set_size(nx, ny, nz);
		sigma_image = params.set_default("sigma", (EMData*)0);
		ignore0 = params["ignore0"];

		normimage = params.set_default("normimage", (EMData*)0);
		if (ignore0 && normimage==0) { normimage=new EMData(nx,ny,nz); freenorm=1; }
		if (normimage) normimage->to_zero();
	}

	float *result_data = result->get_data();
	float *sigma_image_data = 0;
	if (sigma_image) {
		sigma_image->set_size(nx, ny, nz);
		sigma_image_data = sigma_image->get_data();
	}

	float * image_data = image->get_data();

	if (!ignore0) {
		for (size_t j = 0; j < image_size; ++j) {
			float f = image_data[j];
			result_data[j] += f;
			if (sigma_image_data) {
				sigma_image_data[j] += f * f;
			}
		}
	}
	else {
		for (size_t j = 0; j < image_size; ++j) {
			float f = image_data[j];
			if (f) {
				result_data[j] += f;
				if (sigma_image_data) {
					sigma_image_data[j] += f * f;
				}
				normimage->set_value_at_fast(j,normimage->get_value_at(j)+1.0);
			}
		}
	}
}

EMData * ImageAverager::finish()
{
	if (result && nimg > 1) {
		size_t image_size = (size_t)result->get_xsize() * result->get_ysize() * result->get_zsize();
		float * result_data = result->get_data();

		if (!ignore0) {
			for (size_t j = 0; j < image_size; ++j) {
				result_data[j] /= nimg;
			}

			if (sigma_image) {
				float * sigma_image_data = sigma_image->get_data();

				for (size_t j = 0; j < image_size; ++j) {
					float f1 = sigma_image_data[j] / nimg;
					float f2 = result_data[j];
					sigma_image_data[j] = sqrt(f1 - f2 * f2);
				}

				sigma_image->update();
			}
		}
		else {
			for (size_t j = 0; j < image_size; ++j) {
				if (normimage->get_value_at(j)>0) result_data[j] /= normimage->get_value_at(j);
			}
			if (sigma_image) {
				float * sigma_image_data = sigma_image->get_data();

				for (size_t j = 0; j < image_size; ++j) {
					float f1 = 0;
					if (normimage->get_value_at(j)>0) f1=sigma_image_data[j] / normimage->get_value_at(j);
					float f2 = result_data[j];
					sigma_image_data[j] = sqrt(f1 - f2 * f2);
				}

				sigma_image->update();
			}
		}

		result->update();

	}		
	result->set_attr("ptcl_repr",nimg);

	if (freenorm) { delete normimage; normimage=(EMData*)0; }

	return result;
}

VarianceAverager::VarianceAverager()
	: mean(0), nimg(0)
{

}

void VarianceAverager::add_image(EMData * image)
{
	if (!image) {
		return;
	}

	if (nimg >= 1 && !EMUtil::is_same_size(image, result)) {
		LOGERR("%sAverager can only process same-size Image", get_name().c_str());
		return;
	}

	nimg++;

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	size_t image_size = (size_t)nx * ny * nz;

	if (nimg == 1) {
		result = image->copy_head();
		result->set_size(nx, ny, nz);
	}

	float *result_data = result->get_data();
	float *image_data = image->get_data();
	float *mean_data = mean->get_data();

	for (size_t j = 0; j < image_size; ++j) {
		float f = image_data[j];
		if (f) {
			float delta = f - mean_data[j];
			mean_data[j] += delta / ((float) nimg);
			result_data[j] += delta * (f - mean_data[j]);
		}
	}
}

EMData * VarianceAverager::finish()
{
	if (nimg < 2) {
		LOGERR("Variance calculation requires two or more images");
		return NULL;
	}

	if (result) {
		size_t image_size = (size_t)result->get_xsize() * result->get_ysize() * result->get_zsize();
		float *result_data = result->get_data();
		float tmp = (float)(nimg - 1);
		for (size_t j = 0; j < image_size; ++j) {
			float f = result_data[j];
			if (f) {
				result_data[j] /= tmp;
			}
		}
		result->update();
	}	
	result->set_attr("ptcl_repr",nimg);

	return result;
}

FourierWeightAverager::FourierWeightAverager()
	: normimage(0), freenorm(0), nimg(0)
{

}

void FourierWeightAverager::add_image(EMData * image)
{
	if (!image) {
		return;
	}



	EMData *img=image->do_fft();
	if (nimg >= 1 && !EMUtil::is_same_size(img, result)) {
		LOGERR("%sAverager can only process same-size Image",
			   get_name().c_str());
		return;
	}
	
	nimg++;

	int nx = img->get_xsize();
	int ny = img->get_ysize();
	int nz = 1;
//	size_t image_size = (size_t)nx * ny * nz;

	XYData *weight=(XYData *)image->get_attr("avg_weight");
	
	if (nimg == 1) {
		result = new EMData(nx,ny,nz);
		result->set_complex(true);
		result->to_zero();

		normimage = params.set_default("normimage", (EMData*)0);
		if (normimage==0) { normimage=new EMData(nx/2,ny,nz); freenorm=1; }
		normimage->to_zero();
	}

	// We're using routines that handle complex image wraparound for us, so we iterate over the half-plane
	for (int y=-ny/2; y<ny/2; y++) {
		for (int x=0; x<nx/2; x++) {
			std::complex<float> v=img->get_complex_at(x,y);
			float r=Util::hypot2(y/(float)ny,x/(float)nx);
			float wt=weight->get_yatx(r);
			result->set_complex_at(x,y,result->get_complex_at(x,y)+v*wt);
			normimage->set_value_at(x,y+ny/2,normimage->get_value_at(x,y+ny/2)+wt);
		}
	}

	delete img;
}

EMData * FourierWeightAverager::finish()
{
	EMData *ret = (EMData *)0;
	
	if (result && nimg >= 1) {
	// We're using routines that handle complex image wraparound for us, so we iterate over the half-plane
		int nx = result->get_xsize();
		int ny = result->get_ysize();
		
		for (int y=-ny/2; y<ny/2; y++) {
			for (int x=0; x<nx/2; x++) {
				float norm=normimage->get_value_at(x,y+ny/2);
				if (norm<=0) result->set_complex_at(x,y,0.0f);
				else result->set_complex_at(x,y,result->get_complex_at(x,y)/norm);
			}
		}

		result->update();
//		result->mult(1.0f/(float)result->get_attr("sigma"));
//		result->write_image("tmp.hdf",0);
//		printf("%g %g %g\n",(float)result->get_attr("sigma"),(float)result->get_attr("minimum"),(float)result->get_attr("maximum"));
		ret=result->do_ift();
		delete result;
		result=(EMData*) 0;
	}
	ret->set_attr("ptcl_repr",nimg);

	if (freenorm) { delete normimage; normimage=(EMData*)0; }
	nimg=0;

	return ret;
}

LocalWeightAverager::LocalWeightAverager()
	: normimage(0), freenorm(0), nimg(0)
{

}

void LocalWeightAverager::add_image(EMData * image)
{
	if (!image) {
		return;
	}
	
	nimg++;

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	if (nimg == 1) {
		result = new EMData(nx,ny,nz);
		result->to_zero();

		normimage = params.set_default("normimage", (EMData*)0);
		if (normimage==0) { normimage=new EMData(nx,ny,nz); freenorm=1; }
		normimage->to_zero();
		normimage->add(0.0000001f);
	}

	images.push_back(image->copy());
}

EMData * LocalWeightAverager::finish()
{
	if (nimg==0) return NULL;
	
	int nx = images.front()->get_xsize();
	int ny = images.front()->get_ysize();
	int nz = images.front()->get_zsize();
	
	EMData *ret = new EMData(nx,ny,nz);
	EMData *stg1 = new EMData(nx,ny,nz);
	
	for (std::vector<EMData*>::iterator im = images.begin(); im!=images.end(); ++im) stg1->add(**im);
	stg1->process_inplace("normalize.edgemean");
	
//	std::vector<EMData*> weights;
	for (std::vector<EMData*>::iterator im = images.begin(); im!=images.end(); ++im) {
		EMData *imc=(*im)->copy();
		imc->mult(*stg1);
		imc->process_inplace("filter.lowpass.gauss",Dict("cutoff_freq",0.02f));
		imc->process_inplace("threshold.belowtozero",Dict("minval",0.0f));
//		imc->process_inplace("math.sqrt");
//		imc->process_inplace("math.pow",Dict("pow",0.25));
		(*im)->mult(*imc);
		result->add(**im);
		normimage->add(*imc);
		delete *im;
	}
	
	result->div(*normimage);
	
	result->set_attr("ptcl_repr",nimg);
	
	if (freenorm) { delete normimage; normimage=(EMData*)0; }
	nimg=0;

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
	size_t image_size = (size_t)nx * ny * nz;

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
		for (size_t j = 0; j < image_size; ++j) {
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

			for (size_t j = 0; j < image_size; ++j) {
				sigma_image_data[j] *= sigma_image_data[j];
			}
		}

		image0->update();
		memcpy(result_data, image0_data, image_size * sizeof(float));

		for (size_t i = 1; i < image_list.size(); i++) {
			EMData *image = image_list[i];

			if (EMUtil::is_same_size(image, result)) {
				float *image_data = image->get_data();

				for (size_t j = 0; j < image_size; ++j) {
					result_data[j] += image_data[j];
				}

				if (sigma_image_data) {
					for (size_t j = 0; j < image_size; ++j) {
						sigma_image_data[j] += image_data[j] * image_data[j];
					}
				}

				image->update();
				c++;
			}
		}

		for (size_t j = 0; j < image_size; ++j) {
			result_data[j] /= static_cast < float >(c);
		}
	}

	if (sigma_image_data) {
		for (size_t j = 0; j < image_size; ++j) {
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
	EMData *owner = params.set_default("owner", (EMData*)0);

	float thisown = image->get_attr_default("ortid",(float)nimg);
	nimg++;

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	if (nimg == 1) {
		result = image->copy();
		max = params["max"];
		return;
	}

	if (max) {
		for (int z=0; z<nz; z++) {
			for (int y=0; y<ny; y++) {
				for (int x=0; x<nx; x++) {
					if (result->get_value_at(x,y,z)<=image->get_value_at(x,y,z)) {
						result->set_value_at(x,y,z,image->get_value_at(x,y,z));
						if (owner) owner->set_value_at(x,y,z,thisown);
					}
				}
			}
		}
	}
	else {
		for (int z=0; z<nz; z++) {
			for (int y=0; y<ny; y++) {
				for (int x=0; x<nx; x++) {
					if (result->get_value_at(x,y,z)>image->get_value_at(x,y,z)) {
						result->set_value_at(x,y,z,image->get_value_at(x,y,z)); 
						if (owner) owner->set_value_at(x,y,z,thisown);
					}
				}
			}
		}
	}

}

EMData *MinMaxAverager::finish()
{
	result->update();
	result->set_attr("ptcl_repr",nimg);
	
	if (result && nimg > 1) return result;

	return NULL;
}

AbsMaxMinAverager::AbsMaxMinAverager() : nimg(0)
{
	/*move max out of initializer list, since this max(0) is considered as a macro
	 * in Visual Studio, which we defined somewhere else*/
	min = 0;
}

void AbsMaxMinAverager::add_image(EMData * image)
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

	size_t imgsize = (size_t)nx*ny*nz;

	if (nimg == 1) {
		result = image->copy();
		min = params["min"];
		return;
	}

	float * data 	 = result->get_data();
	float * src_data = image->get_data();

	for(size_t i=0; i<imgsize; ++i) {
		if(!min) {	//average to maximum by default
			if (fabs(data[i]) < fabs(src_data[i])) data[i]=src_data[i];
		}
		else {	//average to minimum if set 'min'
			if (fabs(data[i]) > fabs(src_data[i])) data[i]=src_data[i];
		}
	}
}

EMData *AbsMaxMinAverager::finish()
{
	result->update();
	result->set_attr("ptcl_repr",nimg);

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
	size_t image_size = (size_t)nx * ny * nz;

	if (nimg == 1) {
		result = image->copy_head();
		result->set_size(nx, ny, nz);
		sigma_image = image->copy_head();
		sigma_image->set_size(nx, ny, nz);
	}

	float *image_data = image->get_data();
	float *result_data = result->get_data();
	float *sigma_image_data = sigma_image->get_data();

	for (size_t j = 0; j < image_size; ++j) {
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
	size_t image_size = (size_t)nx * ny * nz;

	float *result_data = result->get_data();
	float *sigma_image_data = sigma_image->get_data();

	for (size_t j = 0; j < image_size; ++j) {
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
//string cc=ctf->to_string();
//FILE *out=fopen("ctf.txt","a");
//fprintf(out,"%s\n",cc.c_str());
//fclose(out);
	float b=ctf->bfactor;
	ctf->bfactor=100.0;			// FIXME - this is a temporary fixed B-factor which does a (very) little sharpening

//	if (nimg==1) unlink("snr.hdf");

	EMData *snr = result -> copy();
	ctf->compute_2d_complex(snr,Ctf::CTF_SNR);
//	snr->write_image("snr.hdf",-1);
	EMData *ctfi = result-> copy();
	ctf->compute_2d_complex(ctfi,Ctf::CTF_AMP);

	ctf->bfactor=b;	// return to its original value

	float *outd = result->get_data();
	float *ind = fft->get_data();
	float *snrd = snr->get_data();
	float *ctfd = ctfi->get_data();

	size_t sz=snr->get_xsize()*snr->get_ysize();
	for (size_t i = 0; i < sz; i+=2) {
		if (snrd[i]<0) snrd[i]=0.001;	// Used to be 0. See ctfcauto averager
		ctfd[i]=fabs(ctfd[i]);
		if (ctfd[i]<.05) ctfd[i]=0.05f;
//		{
//			if (snrd[i]<=0) ctfd[i]=.05f;
//			else ctfd[i]=snrd[i]*10.0f;
//		}
		outd[i]+=ind[i]*snrd[i]/ctfd[i];
		outd[i+1]+=ind[i+1]*snrd[i]/ctfd[i];
	}

	if (nimg==1) {
		snrsum=snr->copy_head();
		float *ssnrd=snrsum->get_data();
		// we're only using the real component, and we need to start with 1.0
		for (size_t i = 0; i < sz; i+=2) { ssnrd[i]=1.0; ssnrd[i+1]=0.0; }
	}
//	snr->write_image("snr.hdf",-1);
	snr->process_inplace("math.absvalue");
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

//snrsum->write_image("snrsum.hdf",-1);
	//size_t sz=result->get_xsize()*result->get_ysize();
	int nx=result->get_xsize();
	int ny=result->get_ysize();	
	float *snrsd=snrsum->get_data();
	float *outd=result->get_data();

	int rm=(ny-2)*(ny-2)/4;
	for (int j=0; j<ny; j++) {
		for (int i=0; i<nx; i+=2) {
			size_t ii=i+j*nx;
			if ((j<ny/2 && i*i/4+j*j>rm) ||(j>=ny/2 && i*i/4+(ny-j)*(ny-j)>rm) || snrsd[ii]==0) { outd[ii]=outd[ii+1]=0; continue; }
			outd[ii]/=snrsd[ii];		// snrsd contains total SNR
			outd[ii+1]/=snrsd[ii];
		}
	}

	result->update();
	result->set_attr("ptcl_repr",nimg);
	result->set_attr("ctf_snr_total",snrsum->calc_radial_dist(snrsum->get_ysize()/2,0,1,false));
	result->set_attr("ctf_wiener_filtered",1);
	
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
		if (snrd[i]<=0) snrd[i]=0.001;		// used to be 0. Trying to insure that there is always at least a little signal used. In cases with dense particles, SNR may be dramatically underestimated
		ctfd[i]=fabs(ctfd[i]);
		
		// This limits the maximum possible amplification in CTF correction to 10x
		if (ctfd[i]<.05)  ctfd[i]=0.05f;
//		{
//			if (snrd[i]<=0) ctfd[i]=.05f;
//			else ctfd[i]=snrd[i]*10.0f;
//		}
		
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
	snr->process_inplace("math.absvalue");
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
	//size_t sz=result->get_xsize()*result->get_ysize();
	int nx=result->get_xsize();
	int ny=result->get_ysize();	
	float *snrsd=snrsum->get_data();
	float *outd=result->get_data();

	int rm=(ny-2)*(ny-2)/4;
	for (int j=0; j<ny; j++) {
		for (int i=0; i<nx; i+=2) {
			size_t ii=i+j*nx;
			if ((j<ny/2 && i*i/4+j*j>rm) ||(j>=ny/2 && i*i/4+(ny-j)*(ny-j)>rm) || snrsd[ii]==0) { outd[ii]=outd[ii+1]=0; continue; }
			// we aren't wiener filtering, but if the total SNR is too low, we don't want TOO much exaggeration of noise
			if (snrsd[ii]<.05) {		
				outd[ii]*=20.0;		// 1/0.05
				outd[ii+1]*=20.0;
			}
			else {
				outd[ii]/=snrsd[ii];		// snrsd contains total SNR
				outd[ii+1]/=snrsd[ii];
			}
		}
	}
	result->update();
	result->set_attr("ptcl_repr",nimg);
	result->set_attr("ctf_snr_total",snrsum->calc_radial_dist(snrsum->get_ysize()/2,0,1,false));
	result->set_attr("ctf_wiener_filtered",0);
	
/*	snrsum->write_image("snr.hdf",-1);
	result->write_image("avg.hdf",-1);*/
	
	delete snrsum;
	EMData *ret=result->do_ift();
	delete result;
	result=NULL;
	return ret;
}

CtfWtAverager::CtfWtAverager()
	: nimg(0)
{

}


void CtfWtAverager::add_image(EMData * image)
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

	EMData *ctfi = result-> copy();
	float b=ctf->bfactor;
	ctf->bfactor=0;		// no B-factor used in weight
	ctf->compute_2d_complex(ctfi,Ctf::CTF_INTEN);
	ctf->bfactor=b;	// return to its original value

	float *outd = result->get_data();
	float *ind = fft->get_data();
	float *ctfd = ctfi->get_data();

	size_t sz=ctfi->get_xsize()*ctfi->get_ysize();
	for (size_t i = 0; i < sz; i+=2) {
		
		// CTF weight
		outd[i]+=ind[i]*ctfd[i];
		outd[i+1]+=ind[i+1]*ctfd[i];
	}

	if (nimg==1) {
		ctfsum=ctfi->copy_head();
		ctfsum->to_zero();
	}
	ctfsum->add(*ctfi);

	delete ctf;
	delete fft;
	delete ctfi;
}

EMData * CtfWtAverager::finish()
{
	int nx=result->get_xsize();
	int ny=result->get_ysize();	
	float *ctfsd=ctfsum->get_data();
	float *outd=result->get_data();

	for (int j=0; j<ny; j++) {
		for (int i=0; i<nx; i+=2) {
			size_t ii=i+j*nx;
			outd[ii]/=ctfsd[ii];		// snrsd contains total SNR
			outd[ii+1]/=ctfsd[ii];
		}
	}
	result->update();
	result->set_attr("ptcl_repr",nimg);
//	result->set_attr("ctf_total",ctfsum->calc_radial_dist(ctfsum->get_ysize()/2,0,1,false));
	result->set_attr("ctf_wiener_filtered",0);
	
/*	snrsum->write_image("snr.hdf",-1);
	result->write_image("avg.hdf",-1);*/
	
	delete ctfsum;
	EMData *ret=result->do_ift();
	delete result;
	result=NULL;
	return ret;
}

CtfWtFiltAverager::CtfWtFiltAverager()
{
	nimg[0]=0;
	nimg[1]=0;
	eo=-1;
}


void CtfWtFiltAverager::add_image(EMData * image)
{
	if (!image) {
		return;
	}



	EMData *fft=image->do_fft();

	if (nimg[0] >= 1 && !EMUtil::is_same_size(fft, results[0])) {
		LOGERR("%s Averager can only process images of the same size", get_name().c_str());
		return;
	}

	if (eo==-1) {
		results[0] = fft->copy_head();
		results[0]->to_zero();
		results[1] = fft->copy_head();
		results[1]->to_zero();
		eo=1;
	}

	eo^=1;
	nimg[eo]++;

	
	EMData *ctfi = results[0]-> copy();
	if (image->has_attr("ctf")) {
		Ctf *ctf = (Ctf *)image->get_attr("ctf");

		float b=ctf->bfactor;
		ctf->bfactor=0;		// no B-factor used in weight, not strictly threadsafe, but shouldn't be a problem
		ctf->compute_2d_complex(ctfi,Ctf::CTF_INTEN);
		ctf->bfactor=b;	// return to its original value
		delete ctf;
	}
	else {
		ctfi->to_one();
	}
		
	float *outd = results[eo]->get_data();
	float *ind = fft->get_data();
	float *ctfd = ctfi->get_data();

	size_t sz=ctfi->get_xsize()*ctfi->get_ysize();
	for (size_t i = 0; i < sz; i+=2) {
		
		// CTF weight
		outd[i]+=ind[i]*ctfd[i];
		outd[i+1]+=ind[i+1]*ctfd[i];
	}

	if (nimg[eo]==1) {
		ctfsum[eo]=ctfi->copy_head();
		ctfsum[eo]->to_zero();
		ctfsum[eo]->add(0.1);		// we start with a value of 0.1 rather than zero to empirically help with situations where the data is incomplete
	}
	ctfsum[eo]->add(*ctfi);

	delete fft;
	delete ctfi;
}

EMData * CtfWtFiltAverager::finish()
{
	if (nimg[0]==0 && nimg[1]==0) return NULL;	// no images
	// Only a single image, so we just return it. No way to filter
	if (nimg[1]==0) {
		EMData *ret=results[0]->do_ift();
		delete results[0];
		delete results[1];
		delete ctfsum[0];
		return ret;
	}

	int nx=results[0]->get_xsize();
	int ny=results[0]->get_ysize();

	for (int k=0; k<2; k++) {
		float *outd=results[k]->get_data();
		float *ctfsd=ctfsum[k]->get_data();
		for (int j=0; j<ny; j++) {
			for (int i=0; i<nx; i+=2) {
				size_t ii=i+j*nx;
				outd[ii]/=ctfsd[ii];
				outd[ii+1]/=ctfsd[ii];
			}
		}
		results[k]->update();
	//	result->set_attr("ctf_total",ctfsum->calc_radial_dist(ctfsum->get_ysize()/2,0,1,false));
		results[0]->set_attr("ctf_wiener_filtered",1);
	}
	
	// compute the Wiener filter from the FSC
	std::vector<float> fsc=results[0]->calc_fourier_shell_correlation(results[1]);
	int third=fsc.size()/3;
	for (int i=third; i<third*2; i++) {
		if (fsc[i]>=.9999) fsc[i]=.9999;
		if (fsc[i]<.001) fsc[i]=.001;
		float snr=fsc[i]/(1.0-fsc[i]);
		fsc[i]=snr*snr/(snr*snr+1.0);
	}
	
	
	results[0]->add(*results[1]);
	
	float c;
	for (int j=-ny/2; j<ny/2; j++) {
		for (int i=0; i<nx/2; i++) {
			int r=(int)Util::hypot_fast(i,j);
			if (r>=third) c=0.0;
			else c=fsc[third+r];
			results[0]->set_complex_at(i,j,results[0]->get_complex_at(i,j)*c);
		}
	}
	
	EMData *ret=results[0]->do_ift();
	ret->set_attr("ptcl_repr",nimg[0]+nimg[1]);
	
/*	snrsum->write_image("snr.hdf",-1);
	result->write_image("avg.hdf",-1);*/
	
	delete ctfsum[0];
	delete ctfsum[1];
	delete results[0];
	delete results[1];
	results[0]=results[1]=NULL;
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
	size_t image_size = (size_t)nx * ny * nz;

	EMData *result = new EMData();
	result->set_size(nx, ny, nz);

	EMData *sigma_image = new EMData();
	sigma_image->set_size(nx, ny, nz);

	float *image0_data = image0->get_data();
	float *result_data = result->get_data();
	float *sigma_image_data = sigma_image->get_data();

	memcpy(result_data, image0_data, image_size * sizeof(float));
	memcpy(sigma_image_data, image0_data, image_size * sizeof(float));

	for (size_t j = 0; j < image_size; ++j) {
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
	for (size_t j = 0; j < image_size; ++j) {
		float f1 = sigma_image_data[j] / c;
		float f2 = result_data[j] / c;
		sigma_image_data[j] = sqrt(f1 - f2 * f2) / sqrt(c);
	}


	for (size_t j = 0; j < image_size; ++j) {
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




void EMAN::dump_averagers()
{
	dump_factory < Averager > ();
}

map<string, vector<string> > EMAN::dump_averagers_list()
{
	return dump_factory_list < Averager > ();
}
