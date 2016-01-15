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
#include "emfft.h"
#include "cmp.h"
#include "aligner.h"
#include "averager.h"
#include "emdata.h"
#include "processor.h"
#include "util.h"
#include "symmetry.h"
#include <gsl/gsl_multimin.h>
#include "plugins/aligner_template.h"

#ifdef EMAN2_USING_CUDA
	#include "cuda/cuda_processor.h"
	#include "cuda/cuda_cmp.h"
#endif

#ifdef SPARX_USING_CUDA
	#include <sparx/cuda/cuda_ccf.h>
#endif

#define EMAN2_ALIGNER_DEBUG 0

using namespace EMAN;

const string TranslationalAligner::NAME = "translational";
const string RotationalAligner::NAME = "rotational";
const string RotationalAlignerIterative::NAME = "rotational_iterative";
const string RotatePrecenterAligner::NAME = "rotate_precenter";
const string RotateTranslateAligner::NAME = "rotate_translate";
const string RotateTranslateScaleAligner::NAME = "rotate_translate_scale";
const string RotateTranslateAlignerIterative::NAME = "rotate_translate_iterative";
const string RotateTranslateScaleAlignerIterative::NAME = "rotate_trans_scale_iter";
const string RotateTranslateAlignerPawel::NAME = "rotate_translate_resample";
const string RotateTranslateBestAligner::NAME = "rotate_translate_best";
const string RotateFlipAligner::NAME = "rotate_flip";
const string RotateFlipAlignerIterative::NAME = "rotate_flip_iterative";
const string RotateTranslateFlipAligner::NAME = "rotate_translate_flip";
const string RotateTranslateFlipScaleAligner::NAME = "rotate_trans_flip_scale";
const string RotateTranslateFlipAlignerIterative::NAME = "rotate_translate_flip_iterative";
const string RotateTranslateFlipScaleAlignerIterative::NAME = "rotate_trans_flip_scale_iter";
const string RotateTranslateFlipAlignerPawel::NAME = "rotate_translate_flip_resample";
const string RTFExhaustiveAligner::NAME = "rtf_exhaustive";
const string RTFSlowExhaustiveAligner::NAME = "rtf_slow_exhaustive";
const string RefineAligner::NAME = "refine";
const string RefineAlignerCG::NAME = "refinecg";
const string SymAlignProcessorQuat::NAME = "symalignquat";
const string SymAlignProcessor::NAME = "symalign";
const string Refine3DAlignerGrid::NAME = "refine_3d_grid";
const string Refine3DAlignerQuaternion::NAME = "refine_3d";
const string RT3DGridAligner::NAME = "rotate_translate_3d_grid";
const string RT3DSphereAligner::NAME = "rotate_translate_3d";
const string RT3DTreeAligner::NAME = "rotate_translate_3d_tree";
const string RT3DSymmetryAligner::NAME = "rotate_symmetry_3d";
const string FRM2DAligner::NAME = "frm2d";
const string ScaleAligner::NAME = "scale";


template <> Factory < Aligner >::Factory()
{
	force_add<TranslationalAligner>();
	force_add<RotationalAligner>();
	force_add<RotationalAlignerIterative>();
	force_add<RotatePrecenterAligner>();
	force_add<RotateTranslateAligner>();
	force_add<RotateTranslateScaleAligner>();
	force_add<RotateTranslateAlignerIterative>();
	force_add<RotateTranslateScaleAlignerIterative>();
	force_add<RotateTranslateAlignerPawel>();
	force_add<RotateFlipAligner>();
	force_add<RotateFlipAlignerIterative>();
	force_add<RotateTranslateFlipAligner>();
	force_add<RotateTranslateFlipScaleAligner>();
	force_add<RotateTranslateFlipAlignerIterative>();
	force_add<RotateTranslateFlipScaleAlignerIterative>();
	force_add<RotateTranslateFlipAlignerPawel>();
	force_add<RTFExhaustiveAligner>();
	force_add<RTFSlowExhaustiveAligner>();
	force_add<SymAlignProcessor>();
	force_add<RefineAligner>();
	force_add<RefineAlignerCG>();
	force_add<SymAlignProcessorQuat>();
	force_add<Refine3DAlignerGrid>();
	force_add<Refine3DAlignerQuaternion>();
	force_add<RT3DGridAligner>();
	force_add<RT3DSphereAligner>();
	force_add<RT3DTreeAligner>();
	force_add<RT3DSymmetryAligner>();
	force_add<FRM2DAligner>();
	force_add<ScaleAligner>();
//	force_add<XYZAligner>();
}

vector<Dict> Aligner::xform_align_nbest(EMData *, EMData *, const unsigned int, const string &, const Dict&) const
{
	vector<Dict> solns;
	return solns;
}

EMData* ScaleAlignerABS::align_using_base(EMData * this_img, EMData * to,
			const string & cmp_name, const Dict& cmp_params) const
{
	//get the scale range
	float min =  params.set_default("min",0.95f);
	float max = params.set_default("max",1.05f);
	float step = params.set_default("step",0.01f);

	// crate the starting transform
	Transform t = Transform();
	t.set_scale(max);

	//save orignal data
	float* oridata = this_img->get_data();

	//get the transform processor and cast to correct factory product
	Processor* proc = Factory <Processor>::get("xform", Dict());
	TransformProcessor* xform = dynamic_cast<TransformProcessor*>(proc);

	// Using the following method we only create one EMdata object. If I just used the processor, then I would create many EMdata objects
	EMData* result = 0;
//	float bestscore = numeric_limits<float>::infinity();
	float bestscore = 1.0e37;

	for(float i = max; i > min; i-=step){

		//scale the image
		float* des_data = xform->transform(this_img,t);
		this_img->set_data(des_data);
		this_img->update();

		//check compairsion
		EMData* aligned = this_img->align(basealigner, to, basealigner_params, cmp_name, cmp_params);
		float score = aligned->cmp(cmp_name, to, cmp_params);
		if(score < bestscore){
			bestscore = score;
			//If we just reassign w/o cleaing up we will get memory leaks!
			if(result != 0) delete result;
			result = aligned;
			result->set_attr("scalefactor",i);
		}else{
			delete aligned;
		}
		//clean up scaled image data
		delete des_data;

		t.set_scale(i);

		//reset original data
		this_img->set_data(oridata);
	}

	if (!result) throw UnexpectedBehaviorException("Alignment score is infinity! Something is seriously wrong with the data!");
	if (proc != 0) delete proc;

	return result;

};

EMData* ScaleAligner::align(EMData * this_img, EMData *to,
			const string& cmp_name, const Dict& cmp_params) const
{

	//get the scale range
	float min =  params.set_default("min",0.95f);
	float max = params.set_default("max",1.05f);
	float step = params.set_default("step",0.01f);

	Transform t = Transform();
	t.set_scale(max);
	float* oridata = this_img->get_data();

	//get the transform processor and cast to correct factory product
	Processor* proc = Factory <Processor>::get("xform", Dict());
	TransformProcessor* xform = dynamic_cast<TransformProcessor*>(proc);

	// Using the following method we only create one EMdata object. If I just used the processor, then I would create many EMdata objects
	float bestscale = 1.0;
	float bestscore = 1.0e37;

	for(float i = max; i > min; i-=step){

		float* des_data = xform->transform(this_img,t);
		this_img->set_data(des_data);
		this_img->update();

		//check compairsion
		float score = this_img->cmp(cmp_name, to, cmp_params);
		if(score < bestscore){
			bestscore = score;
			bestscale = i;
		}
		//clean up scaled image data
		delete des_data;

		t.set_scale(i);

		//reset original data
		this_img->set_data(oridata);
	}



	//Return scaled image
	t.set_scale(bestscale);
	EMData* result = this_img->process("xform",Dict("transform",&t));
	result->set_attr("scalefactor",bestscale);
	if (proc != 0) delete proc;

	return result;

}

// Note, the translational aligner assumes that the correlation image
// generated by the calc_ccf function is centered on the bottom left corner
// That is, if you did at calc_cff using identical images, the
// peak would be at 0,0
EMData *TranslationalAligner::align(EMData * this_img, EMData *to,
					const string&, const Dict&) const
{
	if (!this_img) {
		return 0;
	}

	if (to && !EMUtil::is_same_size(this_img, to))
		throw ImageDimensionException("Images must be the same size to perform translational alignment");

	EMData *cf = 0;
	int nx = this_img->get_xsize();
	int ny = this_img->get_ysize();
	int nz = this_img->get_zsize();

	int masked = params.set_default("masked",0);
	int useflcf = params.set_default("useflcf",0);
	bool use_cpu = true;

#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		//if(!this_img->getcudarwdata()) this_img->copy_to_cuda();
		//if(to && !to->getcudarwdata()) to->copy_to_cuda();
		//if (masked) throw UnexpectedBehaviorException("Masked is not yet supported in CUDA");
		//if (useflcf) throw UnexpectedBehaviorException("Useflcf is not yet supported in CUDA");
 		//cout << "Translate on GPU" << endl;
		//use_cpu = false;
		//cf = this_img->calc_ccf(to);
	}
#endif // EMAN2_USING_CUDA

	if (use_cpu) {
		if (useflcf) cf = this_img->calc_flcf(to);
		else cf = this_img->calc_ccf(to);
	}
	//return cf;
	// This is too expensive, esp for CUDA(we we can fix later
	if (masked) {
		EMData *msk=this_img->process("threshold.notzero");
		EMData *sqr=to->process("math.squared");
		EMData *cfn=msk->calc_ccf(sqr);
		cfn->process_inplace("math.sqrt");
		float *d1=cf->get_data();
		float *d2=cfn->get_data();
		for (size_t i=0; i<(size_t)nx*ny*nz; ++i) {
			if (d2[i]!=0) d1[i]/=d2[i];
		}
		cf->update();
		delete msk;
		delete sqr;
		delete cfn;
	}

	int maxshiftx = params.set_default("maxshift",-1);
	int maxshifty = params["maxshift"];
	int maxshiftz = params["maxshift"];
	int nozero = params["nozero"];

	if (maxshiftx <= 0) {
		maxshiftx = nx / 4;
		maxshifty = ny / 4;
		maxshiftz = nz / 4;
	}

	if (maxshiftx > nx / 2 - 1) maxshiftx = nx / 2 - 1;
	if (maxshifty > ny / 2 - 1)	maxshifty = ny / 2 - 1;
	if (maxshiftz > nz / 2 - 1) maxshiftz = nz / 2 - 1;

	if (nx == 1) maxshiftx = 0; // This is justhere for completeness really... plus it saves errors
	if (ny == 1) maxshifty = 0;
	if (nz == 1) maxshiftz = 0;

	// If nozero the portion of the image in the center (and its 8-connected neighborhood) is zeroed
	if (nozero) {
		cf->zero_corner_circulant(1);
	}

	IntPoint peak;
#ifdef EMAN2_USING_CUDA
	if (!use_cpu) {
		cout << "USe CUDA TA 2" << endl;
		if (nozero) throw UnexpectedBehaviorException("Nozero is not yet supported in CUDA");
		CudaPeakInfo* data = calc_max_location_wrap_cuda(cf->getcudarwdata(), cf->get_xsize(), cf->get_ysize(), cf->get_zsize(), maxshiftx, maxshifty, maxshiftz);
		peak = IntPoint(data->px,data->py,data->pz);
		free(data);
	}
#endif // EMAN2_USING_CUDA

	if (use_cpu) {
		peak = cf->calc_max_location_wrap(maxshiftx, maxshifty, maxshiftz);
	}
	//cout << -peak[0] << " " << -peak[1] << " " << -peak[2] << endl;
	Vec3f cur_trans = Vec3f ( (float)-peak[0], (float)-peak[1], (float)-peak[2]);
	//cout << peak[0] << " " << peak[1] << endl;

	if (!to) {
		cur_trans /= 2.0f; // If aligning theimage to itself then only go half way -
		int intonly = params.set_default("intonly",false);
		if (intonly) {
			cur_trans[0] = floor(cur_trans[0] + 0.5f);
			cur_trans[1] = floor(cur_trans[1] + 0.5f);
			cur_trans[2] = floor(cur_trans[2] + 0.5f);
		}
	}

	if( cf ){
		delete cf;
		cf = 0;
	}

	Dict params("trans",static_cast< vector<int> >(cur_trans));
	if (use_cpu){
		cf=this_img->process("xform.translate.int",params);
	}
	Transform t;
	t.set_trans(cur_trans);

#ifdef EMAN2_USING_CUDA
	if (!use_cpu) {
		cout << "USe CUDA TA 3" << endl;
		//this will work just fine....
		cf = this_img->process("xform",Dict("transform",&t));
	}
#endif // EMAN2_USING_CUDA

	if ( nz != 1 ) {
//		Transform* t = get_set_align_attr("xform.align3d",cf,this_img);
//		t->set_trans(cur_trans);
		cf->set_attr("xform.align3d",&t);
	} else if ( ny != 1 ) {
		//Transform* t = get_set_align_attr("xform.align2d",cf,this_img);
		cur_trans[2] = 0; // just make sure of it
		t.set_trans(cur_trans);
		cf->set_attr("xform.align2d",&t);
	}
	return cf;
}

EMData * RotationalAligner::align_180_ambiguous(EMData * this_img, EMData * to, int rfp_mode,int zscore) {

	// Make translationally invariant rotational footprints
	EMData* this_img_rfp, * to_rfp;
	if (rfp_mode == 0) {
		this_img_rfp = this_img->make_rotational_footprint_e1();
		to_rfp = to->make_rotational_footprint_e1();
	} else if (rfp_mode == 1) {
		this_img_rfp = this_img->make_rotational_footprint();
		to_rfp = to->make_rotational_footprint();
	} else if (rfp_mode == 2) {
		this_img_rfp = this_img->make_rotational_footprint_cmc();
		to_rfp = to->make_rotational_footprint_cmc();
	} else {
		throw InvalidParameterException("rfp_mode must be 0,1 or 2");
	}
	int this_img_rfp_nx = this_img_rfp->get_xsize();

	// Do row-wise correlation, returning a sum.
	EMData *cf = this_img_rfp->calc_ccfx(to_rfp, 0, this_img->get_ysize(),false,false,zscore);
// cf->process_inplace("normalize");
// cf->write_image("ralisum.hdf",-1);
//
// EMData *cf2 = this_img_rfp->calc_ccfx(to_rfp, 0, this_img->get_ysize(),true);
// cf2->write_image("ralistack.hdf",-1);
// delete cf2;

	// Delete them, they're no longer needed
	delete this_img_rfp; this_img_rfp = 0;
	delete to_rfp; to_rfp = 0;

	// Now solve the rotational alignment by finding the max in the column sum
	float *data = cf->get_data();

	float peak = 0;
	int peak_index = 0;
	Util::find_max(data, this_img_rfp_nx, &peak, &peak_index);

	if( cf ) {
		delete cf;
		cf = 0;
	}
	float rot_angle = (float) (peak_index * 180.0f / this_img_rfp_nx);

	// Return the result
	Transform tmp(Dict("type","2d","alpha",rot_angle));
	cf=this_img->process("xform",Dict("transform",(Transform*)&tmp));
//	Transform* t = get_set_align_attr("xform.align2d",cf,this_img);
//	Dict d("type","2d","alpha",rot_angle);
//	t->set_rotation(d);
	cf->set_attr("xform.align2d",&tmp);
	return cf;
}

EMData *RotationalAligner::align(EMData * this_img, EMData *to,
			const string& cmp_name, const Dict& cmp_params) const
{
	if (!to) throw InvalidParameterException("Can not rotational align - the image to align to is NULL");

#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		//if(!this_img->getcudarwdata()) this_img->copy_to_cuda();
		//if(!to->getcudarwdata()) to->copy_to_cuda();
	}
#endif

	// Perform 180 ambiguous alignment
	int rfp_mode = params.set_default("rfp_mode",2);
	int zscore = params.set_default("zscore",0);
	int ambig180 = params.set_default("ambig180",0);

	EMData* rot_aligned = RotationalAligner::align_180_ambiguous(this_img,to,rfp_mode,zscore);
	Transform * tmp = rot_aligned->get_attr("xform.align2d");
	Dict rot = tmp->get_rotation("2d");
	float rotate_angle_solution = rot["alpha"];
	delete tmp;

	// Don't resolve the 180 degree ambiguity here
	if (ambig180) {
		return rot_aligned;
	}

	EMData *rot_align_180 = rot_aligned->process("math.rotate.180");

	// Generate the comparison metrics for both rotational candidates
	float rot_cmp = rot_aligned->cmp(cmp_name, to, cmp_params);
	float rot_180_cmp = rot_align_180->cmp(cmp_name, to, cmp_params);

	// Decide on the result
	float score = 0.0;
	EMData* result = NULL;
	if (rot_cmp < rot_180_cmp){
		result = rot_aligned;
		score = rot_cmp;
		delete rot_align_180; rot_align_180 = 0;
	} else {
		result = rot_align_180;
		score = rot_180_cmp;
		delete rot_aligned; rot_aligned = 0;
		rotate_angle_solution = rotate_angle_solution-180.0f;
	}

//	Transform* t = get_align_attr("xform.align2d",result);
//	t->set_rotation(Dict("type","2d","alpha",rotate_angle_solution));
	Transform tmp2(Dict("type","2d","alpha",rotate_angle_solution));
	result->set_attr("xform.align2d",&tmp2);
	return result;
}


EMData *RotatePrecenterAligner::align(EMData * this_img, EMData *to,
			const string&, const Dict&) const
{
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
	float a = (float) ((1.0f - 1.0f * peak_index / size) * 180. * 2);

	Transform rot;
	rot.set_rotation(Dict("type","2d","alpha",(float)(a*180./M_PI)));
	EMData* rslt = this_img->process("xform",Dict("transform",&rot));
	rslt->set_attr("xform.align2d",&rot);
//
//	Transform* t = get_set_align_attr("xform.align2d",rslt,this_img);
//	t->set_rotation(Dict("type","2d","alpha",-a));
//
//	EMData* result this_img->transform(Dict("type","2d","alpha",(float)(a*180./M_PI)));
//
//	cf->set_attr("xform.align2d",t);
//	delete t;
//	cf->update();

	if( e1 )
	{
		delete e1;
		e1 = 0;
	}

	if( e2 )
	{
		delete e2;
		e2 = 0;
	}

	if (cf) {
		delete cf;
		cf = 0;
	}
	return rslt;
}

EMData *RotationalAlignerIterative::align(EMData * this_img, EMData *to,
			const string &, const Dict&) const
{
	int r1 = params.set_default("r1",-1);
	int r2 = params.set_default("r2",-1);
	//to start lest try the original size image. If needed, we can pad it....
	EMData * to_polar = to->unwrap(r1,r2,-1,0,0,true);
	EMData * this_img_polar = this_img->unwrap(r1,r2,-1,0,0,true);
	int this_img_polar_nx = this_img_polar->get_xsize();

	EMData * cf = this_img_polar->calc_ccfx(to_polar, 0, this_img->get_ysize());

	//take out the garbage
	delete to_polar; to_polar = 0;
	delete this_img_polar; this_img_polar = 0;

	float * data = cf->get_data();
	float peak = 0;
	int peak_index = 0;
	Util::find_max(data, this_img_polar_nx, &peak, &peak_index);

	delete cf; cf = 0;
	float rot_angle = (float) (peak_index * 360.0f / this_img_polar_nx);

	//return the result
	//cout << rot_angle << endl;
	Transform tmp(Dict("type","2d","alpha",rot_angle));
	EMData * rotimg=this_img->process("xform",Dict("transform",(Transform*)&tmp));
	rotimg->set_attr("xform.align2d",&tmp);

	return rotimg;

}

EMData *RotateTranslateAlignerIterative::align(EMData * this_img, EMData *to,
			const string & cmp_name, const Dict& cmp_params) const
{
	int maxiter = params.set_default("maxiter", 3);

	Dict trans_params;
	trans_params["intonly"] = 0;
	trans_params["maxshift"] = params.set_default("maxshift", -1);
	trans_params["useflcf"] = params.set_default("useflcf",0);
	trans_params["nozero"] = params.set_default("nozero",false);

	Dict rot_params;
	rot_params["r1"] = params.set_default("r1", -1);
	rot_params["r2"] = params.set_default("r2", -1);

	Transform t;
	EMData * moving_img = this_img;
	for(int it = 0; it < maxiter; it++){

		// First do a translational alignment
		EMData * trans_align = moving_img->align("translational", to, trans_params, cmp_name, cmp_params);
		Transform * tt = trans_align->get_attr("xform.align2d");
		t = *tt*t;
		delete tt;

		//now do rotation
		EMData * rottrans_align = trans_align->align("rotational_iterative", to, rot_params, cmp_name, cmp_params);
		Transform * rt = rottrans_align->get_attr("xform.align2d");
		t = *rt*t;
		delete trans_align; trans_align = 0;
		delete rottrans_align; rottrans_align = 0;
		delete rt;

		//this minimizes interpolation errors (all images that are futher processed will be interpolated at most twice)
		if(it > 0){delete moving_img;}

		moving_img = this_img->process("xform",Dict("transform",&t));  //iterate
	}

	//write the total transformation; Avoids interpolation erros
	moving_img->set_attr("xform.align2d", &t);

	return moving_img;
}

EMData *RotateTranslateScaleAlignerIterative::align(EMData * this_img, EMData *to,
			const string & cmp_name, const Dict& cmp_params) const
{

	//Basically copy params into rotate_translate
	basealigner_params["maxshift"] = params.set_default("maxshift", -1);
	basealigner_params["r1"] = params.set_default("r1",-1);
	basealigner_params["r2"] = params.set_default("r2",-1);
	basealigner_params["maxiter"] = params.set_default("maxiter",3);
	basealigner_params["nozero"] = params.set_default("nozero",false);
	basealigner_params["useflcf"] = params.set_default("useflcf",0);

	//return the correct results
	return align_using_base(this_img, to, cmp_name, cmp_params);

}

EMData *RotateTranslateAlignerPawel::align(EMData * this_img, EMData *to,
			const string & cmp_name, const Dict& cmp_params) const
{
	if (cmp_name != "dot" && cmp_name != "ccc") throw InvalidParameterException("Resample aligner only works for dot and ccc");

	int maxtx = params.set_default("tx", 0);
	int maxty = params.set_default("ty", 0);
	int r1 = params.set_default("r1",-1);
	int r2 = params.set_default("r2",-1);

	if(this_img->get_xsize()/2 - 1 - r2 - maxtx <= 0 || (r2 == -1 && maxtx > 0)) throw InvalidParameterException("nx/2 - 1 - r2 - tx must be greater than or = 0");
	if(this_img->get_ysize()/2 - 1 - r2 - maxty <= 0 || (r2 == -1 && maxty > 0)) throw InvalidParameterException("ny/2 - 1 - r2 - ty must be greater than or = 0");

//	float best_peak = -numeric_limits<float>::infinity();
	float best_peak = -1.0e37;
	int best_peak_index = 0;
	int best_tx = 0;
	int best_ty = 0;
	int polarxsize = 0;

	for(int x = -maxtx; x <= maxtx; x++){
		for(int y = -maxty; y <= maxty; y++){

			EMData * to_polar = to->unwrap(r1,r2,-1,0,0,true);
			EMData * this_img_polar = this_img->unwrap(r1,r2,-1,x,y,true);
			EMData * cf = this_img_polar->calc_ccfx(to_polar, 0, this_img_polar->get_ysize());

			polarxsize = this_img_polar->get_xsize();

			//take out the garbage
			delete to_polar; to_polar = 0;
			delete this_img_polar; this_img_polar = 0;

			float *data = cf->get_data();
			float peak = 0;
			int peak_index = 0;
			Util::find_max(data, polarxsize, &peak, &peak_index);
			delete cf; cf = 0;

			if(peak > best_peak) {
				best_peak = peak;
				best_peak_index = peak_index;
				best_tx = x;
				best_ty = y;
			}
		}
	}

	float rot_angle = (float) (best_peak_index * 360.0f / polarxsize);

	//return the result
	Transform tmptt(Dict("type","2d","alpha",0,"tx",-best_tx,"ty",-best_ty));
	Transform tmprot(Dict("type","2d","alpha",rot_angle,"tx",0,"ty",0));
	Transform total = tmprot*tmptt;
	EMData* rotimg=this_img->process("xform",Dict("transform",(Transform*)&total));
	rotimg->set_attr("xform.align2d",&total);

	return rotimg;

}

EMData *RotateTranslateAligner::align(EMData * this_img, EMData *to,
			const string & cmp_name, const Dict& cmp_params) const
{

#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		//if(!this_img->getcudarwdata()) this_img->copy_to_cuda();
		//if(!to->getcudarwdata()) to->copy_to_cuda();
	}
#endif

	// Get the 180 degree ambiguously rotationally aligned and its 180 degree rotation counterpart
	int zscore = params.set_default("zscore",0);
	int rfp_mode = params.set_default("rfp_mode",2);
	EMData *rot_align  =  RotationalAligner::align_180_ambiguous(this_img,to,rfp_mode,zscore);
	Transform * tmp = rot_align->get_attr("xform.align2d");
	Dict rot = tmp->get_rotation("2d");
	float rotate_angle_solution = rot["alpha"];
	delete tmp;

	EMData *rot_align_180 = rot_align->process("math.rotate.180");

	Dict trans_params;
	trans_params["intonly"]  = 0;
	trans_params["maxshift"] = params.set_default("maxshift", -1);
	trans_params["useflcf"]=params.set_default("useflcf",0);

	// Do the first case translational alignment
	trans_params["nozero"]   = params.set_default("nozero",false);
	EMData* rot_trans = rot_align->align("translational", to, trans_params, cmp_name, cmp_params);
	if( rot_align ) { // Clean up
		delete rot_align;
		rot_align = 0;
	}

	// Do the second case translational alignment
	EMData*  rot_180_trans = rot_align_180->align("translational", to, trans_params, cmp_name, cmp_params);
	if( rot_align_180 )	{ // Clean up
		delete rot_align_180;
		rot_align_180 = 0;
	}

	// Finally decide on the result
	float cmp1 = rot_trans->cmp(cmp_name, to, cmp_params);
	float cmp2 = rot_180_trans->cmp(cmp_name, to, cmp_params);

	EMData *result = 0;
	if (cmp1 < cmp2) { // All comparators are defined so default return is "smaller is better"
		if( rot_180_trans )	{
			delete rot_180_trans;
			rot_180_trans = 0;
		}
		result = rot_trans;
	}
	else {
		if( rot_trans )	{
			delete rot_trans;
			rot_trans = 0;
		}
		result = rot_180_trans;
		rotate_angle_solution -= 180.f;
	}

	Transform* t = result->get_attr("xform.align2d");
	t->set_rotation(Dict("type","2d","alpha",rotate_angle_solution));
	result->set_attr("xform.align2d",t);
	delete t;

	return result;
}


EMData *RotateTranslateScaleAligner::align(EMData * this_img, EMData *to,
			const string & cmp_name, const Dict& cmp_params) const
{

	//Basically copy params into rotate_translate
	basealigner_params["maxshift"] = params.set_default("maxshift", -1);
	basealigner_params["rfp_mode"] = params.set_default("rfp_mode",2);
	basealigner_params["useflcf"] = params.set_default("useflcf",0);
	basealigner_params["zscore"] = params.set_default("zscore",0);

	//return the correct results
	return align_using_base(this_img, to, cmp_name, cmp_params);

}

EMData* RotateTranslateFlipAligner::align(EMData * this_img, EMData *to,
										  const string & cmp_name, const Dict& cmp_params) const
{
	// Get the non flipped rotational, tranlsationally aligned image
	Dict rt_params("maxshift", params["maxshift"], "rfp_mode", params.set_default("rfp_mode",2),"useflcf",params.set_default("useflcf",0),"zscore",params.set_default("zscore",0));
	EMData *rot_trans_align = this_img->align("rotate_translate",to,rt_params,cmp_name, cmp_params);

	// Do the same alignment, but using the flipped version of the image
	EMData *flipped = params.set_default("flip", (EMData *) 0);
	bool delete_flag = false;
	if (flipped == 0) {
		flipped = to->process("xform.flip", Dict("axis", "x"));
		delete_flag = true;
	}

	EMData * rot_trans_align_flip = this_img->align("rotate_translate", flipped, rt_params, cmp_name, cmp_params);
	Transform * t = rot_trans_align_flip->get_attr("xform.align2d");
	t->set_mirror(true);
	rot_trans_align_flip->set_attr("xform.align2d",t);
	delete t;

	// Now finally decide on what is the best answer
	float cmp1 = rot_trans_align->cmp(cmp_name, to, cmp_params);
	float cmp2 = rot_trans_align_flip->cmp(cmp_name, flipped, cmp_params);

	if (delete_flag){
		if(flipped) {
			delete flipped;
			flipped = 0;
		}
	}

	EMData *result = 0;
	if (cmp1 < cmp2 )  {

		if( rot_trans_align_flip ) {
			delete rot_trans_align_flip;
			rot_trans_align_flip = 0;
		}
		result = rot_trans_align;
	}
	else {
		if( rot_trans_align ) {
			delete rot_trans_align;
			rot_trans_align = 0;
		}
		result = rot_trans_align_flip;
		result->process_inplace("xform.flip",Dict("axis","x"));
	}

	return result;
}

EMData *RotateTranslateFlipScaleAligner::align(EMData * this_img, EMData *to,
			const string & cmp_name, const Dict& cmp_params) const
{

	//Basically copy params into rotate_translate
	basealigner_params["flip"] = params.set_default("flip", (EMData *) 0);
	basealigner_params["maxshift"] = params.set_default("maxshift", -1);
	basealigner_params["rfp_mode"] = params.set_default("rfp_mode",2);
	basealigner_params["useflcf"] = params.set_default("useflcf",0);
	basealigner_params["zscore"] = params.set_default("zscore",0);

	//return the correct results
	return align_using_base(this_img, to, cmp_name, cmp_params);

}

EMData* RotateTranslateFlipAlignerIterative::align(EMData * this_img, EMData *to,
										  const string & cmp_name, const Dict& cmp_params) const
{
	// Get the non flipped rotational, tranlsationally aligned image
	Dict rt_params("maxshift", params["maxshift"],"r1",params.set_default("r1",-1),"r2",params.set_default("r2",-1));
	EMData *rot_trans_align = this_img->align("rotate_translate_iterative",to,rt_params,cmp_name, cmp_params);

	// Do the same alignment, but using the flipped version of the image
	EMData *flipped = params.set_default("flip", (EMData *) 0);
	bool delete_flag = false;
	if (flipped == 0) {
		flipped = to->process("xform.flip", Dict("axis", "x"));
		delete_flag = true;
	}

	EMData * rot_trans_align_flip = this_img->align("rotate_translate_iterative", flipped, rt_params, cmp_name, cmp_params);
	Transform* t = rot_trans_align_flip->get_attr("xform.align2d");
	t->set_mirror(true);
	rot_trans_align_flip->set_attr("xform.align2d",t);
	delete t;

	// Now finally decide on what is the best answer
	float cmp1 = rot_trans_align->cmp(cmp_name, to, cmp_params);
	float cmp2 = rot_trans_align_flip->cmp(cmp_name, flipped, cmp_params);

	if (delete_flag){
		if(flipped) {
			delete flipped;
			flipped = 0;
		}
	}

	EMData *result = 0;
	if (cmp1 < cmp2 )  {

		if( rot_trans_align_flip ) {
			delete rot_trans_align_flip;
			rot_trans_align_flip = 0;
		}
		result = rot_trans_align;
	}
	else {
		if( rot_trans_align ) {
			delete rot_trans_align;
			rot_trans_align = 0;
		}
		result = rot_trans_align_flip;
		result->process_inplace("xform.flip",Dict("axis","x"));
	}

	return result;
}

EMData *RotateTranslateFlipScaleAlignerIterative::align(EMData * this_img, EMData *to,
			const string & cmp_name, const Dict& cmp_params) const
{

	//Basically copy params into rotate_translate
	basealigner_params["flip"] = params.set_default("flip", (EMData *) 0);
	basealigner_params["maxshift"] = params.set_default("maxshift", -1);
	basealigner_params["r1"] = params.set_default("r1",-1);
	basealigner_params["r2"] = params.set_default("r2",-1);
	basealigner_params["maxiter"] = params.set_default("maxiter",3);

	//return the correct results
	return align_using_base(this_img, to, cmp_name, cmp_params);

}

EMData *RotateTranslateFlipAlignerPawel::align(EMData * this_img, EMData *to,
			const string & cmp_name, const Dict& cmp_params) const
{
	if (cmp_name != "dot" && cmp_name != "ccc") throw InvalidParameterException("Resample aligner only works for dot and ccc");

	int maxtx = params.set_default("tx", 0);
	int maxty = params.set_default("ty", 0);
	int r1 = params.set_default("r1",-1);
	int r2 = params.set_default("r2",-1);

	if(this_img->get_xsize()/2 - 1 - r2 - maxtx <= 0 || (r2 == -1 && maxtx > 0)){
		cout << "\nRunTimeError: nx/2 - 1 - r2 - tx must be greater than or = 0\n" << endl; // For some reason the expection message is not being print, stupid C++
		throw InvalidParameterException("nx/2 - 1 - r2 - tx must be greater than or = 0");
	}
	if(this_img->get_ysize()/2 - 1 - r2 - maxty <= 0 || (r2 == -1 && maxty > 0)){
		cout << "\nRunTimeError:ny/2 - 1 - r2 - ty must be greater than or = 0\n" << endl; // For some reason the expection message is not being print, stupid C++
		throw InvalidParameterException("ny/2 - 1 - r2 - ty must be greater than or = 0");
	}

//	float best_peak = -numeric_limits<float>::infinity();
	float best_peak = -1.0e37;
	int best_peak_index = 0;
	int best_tx = 0;
	int best_ty = 0;
	int polarxsize = 0;
	bool flip = false;

	for(int x = -maxtx; x <= maxtx; x++){
		for(int y = -maxty; y <= maxty; y++){

			EMData * to_polar = to->unwrap(r1,r2,-1,0,0,true);
			EMData * this_img_polar = this_img->unwrap(r1,r2,-1,x,y,true);
			EMData * cfflip = this_img_polar->calc_ccfx(to_polar, 0, this_img_polar->get_ysize(), false, true);
			EMData * cf = this_img_polar->calc_ccfx(to_polar, 0, this_img_polar->get_ysize());

			polarxsize = this_img_polar->get_xsize();

			//take out the garbage
			delete to_polar; to_polar = 0;
			delete this_img_polar; this_img_polar = 0;

			float *data = cf->get_data();
			float peak = 0;
			int peak_index = 0;
			Util::find_max(data, polarxsize, &peak, &peak_index);
			delete cf; cf = 0;

			if(peak > best_peak) {
				best_peak = peak;
				best_peak_index = peak_index;
				best_tx = x;
				best_ty = y;
				flip = false;
			}

			data = cfflip->get_data();
			Util::find_max(data, polarxsize, &peak, &peak_index);
			delete cfflip; cfflip = 0;

			if(peak > best_peak) {
				best_peak = peak;
				best_peak_index = peak_index;
				best_tx = x;
				best_ty = y;
				flip = true;
			}
		}
	}

	float rot_angle = (float) (best_peak_index * 360.0f / polarxsize);

	//return the result
	Transform tmptt(Dict("type","2d","alpha",0,"tx",-best_tx,"ty",-best_ty));
	Transform tmprot(Dict("type","2d","alpha",rot_angle,"tx",0,"ty",0));
	Transform total = tmprot*tmptt;
	EMData * rotimg=this_img->process("xform",Dict("transform",(Transform*)&total));
	rotimg->set_attr("xform.align2d",&total);
	if(flip == true) {
		rotimg->process_inplace("xform.flip",Dict("axis", "x"));
	}

	return rotimg;

}

EMData *RotateFlipAligner::align(EMData * this_img, EMData *to,
			const string& cmp_name, const Dict& cmp_params) const
{
	Dict rot_params("rfp_mode",params.set_default("rfp_mode",2));
	EMData *r1 = this_img->align("rotational", to, rot_params,cmp_name, cmp_params);


	EMData* flipped =to->process("xform.flip", Dict("axis", "x"));
	EMData *r2 = this_img->align("rotational", flipped,rot_params, cmp_name, cmp_params);
	Transform* t = r2->get_attr("xform.align2d");
	t->set_mirror(true);
	r2->set_attr("xform.align2d",t);
	delete t;

	float cmp1 = r1->cmp(cmp_name, to, cmp_params);
	float cmp2 = r2->cmp(cmp_name, flipped, cmp_params);

	delete flipped; flipped = 0;

	EMData *result = 0;

	if (cmp1 < cmp2) {
		if( r2 )
		{
			delete r2;
			r2 = 0;
		}
		result = r1;
	}
	else {
		if( r1 )
		{
			delete r1;
			r1 = 0;
		}
		result = r2;
		result->process_inplace("xform.flip",Dict("axis","x"));
	}

	return result;
}

EMData *RotateFlipAlignerIterative::align(EMData * this_img, EMData *to,
			const string& cmp_name, const Dict& cmp_params) const
{
	Dict rot_params("r1",params.set_default("r1",-1),"r2",params.set_default("r2",-1));
	EMData * r1 = this_img->align("rotational_iterative", to, rot_params,cmp_name, cmp_params);

	EMData * flipped =to->process("xform.flip", Dict("axis", "x"));
	EMData * r2 = this_img->align("rotational_iterative", flipped,rot_params, cmp_name, cmp_params);
	Transform* t = r2->get_attr("xform.align2d");
	t->set_mirror(true);
	r2->set_attr("xform.align2d",t);
	delete t;

	float cmp1 = r1->cmp(cmp_name, to, cmp_params);
	float cmp2 = r2->cmp(cmp_name, flipped, cmp_params);

	delete flipped; flipped = 0;

	EMData *result = 0;

	if (cmp1 < cmp2) {
		if( r2 )
		{
			delete r2;
			r2 = 0;
		}
		result = r1;
	}
	else {
		if( r1 )
		{
			delete r1;
			r1 = 0;
		}
		result = r2;
		result->process_inplace("xform.flip",Dict("axis","x"));
	}

	return result;
}

// David Woolford says FIXME
// You will note the excessive amount of EMData copying that's going in this function
// This is because functions that are operating on the EMData objects are changing them
// and if you do not use copies the whole algorithm breaks. I did not have time to go
// through and rectify this situation.
// David Woolford says - this problem is related to the fact that many functions that
// take EMData pointers as arguments do not take them as constant pointers to constant
// objects, instead they are treated as raw (completely changeable) pointers. This means
// it's hard to track down which functions are changing the EMData objects, because they
// all do (in name). If this behavior is unavoidable then ignore this comment, however if possible it would
// be good to make things const as much as possible. For example in alignment, technically
// the argument EMData objects (raw pointers) should not be altered... should they?
//
// But const can be very annoying sometimes...
EMData *RTFExhaustiveAligner::align(EMData * this_img, EMData *to,
			const string & cmp_name, const Dict& cmp_params) const
{
	EMData *flip = params.set_default("flip", (EMData *) 0);
	int maxshift = params.set_default("maxshift", this_img->get_xsize()/8);
	if (maxshift < 2) throw InvalidParameterException("maxshift must be greater than or equal to 2");

	int ny = this_img->get_ysize();
	int xst = (int) floor(2 * M_PI * ny);
	xst = Util::calc_best_fft_size(xst);

	Dict d("n",2);
	EMData *to_shrunk_unwrapped = to->process("math.medianshrink",d);

	int to_copy_r2 = to_shrunk_unwrapped->get_ysize() / 2 - 2 - maxshift / 2;
	EMData *tmp = to_shrunk_unwrapped->unwrap(4, to_copy_r2, xst / 2, 0, 0, true);
	if( to_shrunk_unwrapped )
	{
		delete to_shrunk_unwrapped;
		to_shrunk_unwrapped = 0;
	}
	to_shrunk_unwrapped = tmp;

	EMData *to_shrunk_unwrapped_copy = to_shrunk_unwrapped->copy();
	EMData* to_unwrapped = to->unwrap(4, to->get_ysize() / 2 - 2 - maxshift, xst, 0, 0, true);
	EMData *to_unwrapped_copy = to_unwrapped->copy();

	bool delete_flipped = true;
	EMData *flipped = 0;
	if (flip) {
		delete_flipped = false;
		flipped = flip;
	}
	else {
		flipped = to->process("xform.flip", Dict("axis", "x"));
	}
	EMData *to_shrunk_flipped_unwrapped = flipped->process("math.medianshrink",d);
	tmp = to_shrunk_flipped_unwrapped->unwrap(4, to_copy_r2, xst / 2, 0, 0, true);
	if( to_shrunk_flipped_unwrapped )
	{
		delete to_shrunk_flipped_unwrapped;
		to_shrunk_flipped_unwrapped = 0;
	}
	to_shrunk_flipped_unwrapped = tmp;
	EMData *to_shrunk_flipped_unwrapped_copy = to_shrunk_flipped_unwrapped->copy();
	EMData* to_flip_unwrapped = flipped->unwrap(4, to->get_ysize() / 2 - 2 - maxshift, xst, 0, 0, true);
	EMData* to_flip_unwrapped_copy = to_flip_unwrapped->copy();

	if (delete_flipped && flipped != 0) {
		delete flipped;
		flipped = 0;
	}

	EMData *this_shrunk_2 = this_img->process("math.medianshrink",d);

	float bestval = FLT_MAX;
	float bestang = 0;
	int bestflip = 0;
	float bestdx = 0;
	float bestdy = 0;

	int half_maxshift = maxshift / 2;

	int ur2 = this_shrunk_2->get_ysize() / 2 - 2 - half_maxshift;
	for (int dy = -half_maxshift; dy <= half_maxshift; dy += 1) {
		for (int dx = -half_maxshift; dx <= half_maxshift; dx += 1) {
#ifdef	_WIN32
			if (_hypot(dx, dy) <= half_maxshift) {
#else
			if (hypot(dx, dy) <= half_maxshift) {
#endif
				EMData *uw = this_shrunk_2->unwrap(4, ur2, xst / 2, dx, dy, true);
				EMData *uwc = uw->copy();
				EMData *a = uw->calc_ccfx(to_shrunk_unwrapped);

				uwc->rotate_x(a->calc_max_index());
				float cm = uwc->cmp(cmp_name, to_shrunk_unwrapped_copy, cmp_params);
				if (cm < bestval) {
					bestval = cm;
					bestang = (float) (2.0 * M_PI * a->calc_max_index() / a->get_xsize());
					bestdx = (float)dx;
					bestdy = (float)dy;
					bestflip = 0;
				}


				if( a )
				{
					delete a;
					a = 0;
				}
				if( uw )
				{
					delete uw;
					uw = 0;
				}
				if( uwc )
				{
					delete uwc;
					uwc = 0;
				}
				uw = this_shrunk_2->unwrap(4, ur2, xst / 2, dx, dy, true);
				uwc = uw->copy();
				a = uw->calc_ccfx(to_shrunk_flipped_unwrapped);

				uwc->rotate_x(a->calc_max_index());
				cm = uwc->cmp(cmp_name, to_shrunk_flipped_unwrapped_copy, cmp_params);
				if (cm < bestval) {
					bestval = cm;
					bestang = (float) (2.0 * M_PI * a->calc_max_index() / a->get_xsize());
					bestdx = (float)dx;
					bestdy = (float)dy;
					bestflip = 1;
				}

				if( a )
				{
					delete a;
					a = 0;
				}

				if( uw )
				{
					delete uw;
					uw = 0;
				}
				if( uwc )
				{
					delete uwc;
					uwc = 0;
				}
			}
		}
	}
	if( this_shrunk_2 )
	{
		delete this_shrunk_2;
		this_shrunk_2 = 0;
	}
	if( to_shrunk_unwrapped )
	{
		delete to_shrunk_unwrapped;
		to_shrunk_unwrapped = 0;
	}
	if( to_shrunk_unwrapped_copy )
	{
		delete to_shrunk_unwrapped_copy;
		to_shrunk_unwrapped_copy = 0;
	}
	if( to_shrunk_flipped_unwrapped )
	{
		delete to_shrunk_flipped_unwrapped;
		to_shrunk_flipped_unwrapped = 0;
	}
	if( to_shrunk_flipped_unwrapped_copy )
	{
		delete to_shrunk_flipped_unwrapped_copy;
		to_shrunk_flipped_unwrapped_copy = 0;
	}
	bestdx *= 2;
	bestdy *= 2;
	bestval = FLT_MAX;

	float bestdx2 = bestdx;
	float bestdy2 = bestdy;
	// Note I tried steps less than 1.0 (sub pixel precision) and it actually appeared detrimental
	// So my advice is to stick with dx += 1.0 etc unless you really are looking to fine tune this
	// algorithm
	for (float dy = bestdy2 - 3; dy <= bestdy2 + 3; dy += 1.0 ) {
		for (float dx = bestdx2 - 3; dx <= bestdx2 + 3; dx += 1.0 ) {

#ifdef	_WIN32
			if (_hypot(dx, dy) <= maxshift) {
#else
			if (hypot(dx, dy) <= maxshift) {
#endif
				EMData *uw = this_img->unwrap(4, this_img->get_ysize() / 2 - 2 - maxshift, xst, (int)dx, (int)dy, true);
				EMData *uwc = uw->copy();
				EMData *a = uw->calc_ccfx(to_unwrapped);

				uwc->rotate_x(a->calc_max_index());
				float cm = uwc->cmp(cmp_name, to_unwrapped_copy, cmp_params);

				if (cm < bestval) {
					bestval = cm;
					bestang = (float)(2.0 * M_PI * a->calc_max_index() / a->get_xsize());
					bestdx = dx;
					bestdy = dy;
					bestflip = 0;
				}

				if( a )
				{
					delete a;
					a = 0;
				}
				if( uw )
				{
					delete uw;
					uw = 0;
				}
				if( uwc )
				{
					delete uwc;
					uwc = 0;
				}
				uw = this_img->unwrap(4, this_img->get_ysize() / 2 - 2 - maxshift, xst, (int)dx, (int)dy, true);
				uwc = uw->copy();
				a = uw->calc_ccfx(to_flip_unwrapped);

				uwc->rotate_x(a->calc_max_index());
				cm = uwc->cmp(cmp_name, to_flip_unwrapped_copy, cmp_params);

				if (cm < bestval) {
					bestval = cm;
					bestang = (float)(2.0 * M_PI * a->calc_max_index() / a->get_xsize());
					bestdx = dx;
					bestdy = dy;
					bestflip = 1;
				}

				if( a )
				{
					delete a;
					a = 0;
				}
				if( uw )
				{
					delete uw;
					uw = 0;
				}
				if( uwc )
				{
					delete uwc;
					uwc = 0;
				}
			}
		}
	}
	if( to_unwrapped ) {delete to_unwrapped;to_unwrapped = 0;}
	if( to_shrunk_unwrapped ) {	delete to_shrunk_unwrapped;	to_shrunk_unwrapped = 0;}
	if (to_unwrapped_copy) { delete to_unwrapped_copy; to_unwrapped_copy = 0; }
	if (to_flip_unwrapped) { delete to_flip_unwrapped; to_flip_unwrapped = 0; }
	if (to_flip_unwrapped_copy) { delete to_flip_unwrapped_copy; to_flip_unwrapped_copy = 0;}

	bestang *= (float)EMConsts::rad2deg;
	Transform t(Dict("type","2d","alpha",(float)bestang));
	t.set_pre_trans(Vec2f(-bestdx,-bestdy));
	if (bestflip) {
		t.set_mirror(true);
	}

	EMData* ret = this_img->process("xform",Dict("transform",&t));
	ret->set_attr("xform.align2d",&t);

	return ret;
}


EMData *RTFSlowExhaustiveAligner::align(EMData * this_img, EMData *to,
			const string & cmp_name, const Dict& cmp_params) const
{

	EMData *flip = params.set_default("flip", (EMData *) 0);
	int maxshift = params.set_default("maxshift", -1);

	EMData *flipped = 0;

	bool delete_flipped = true;
	if (flip) {
		delete_flipped = false;
		flipped = flip;
	}
	else {
		flipped = to->process("xform.flip", Dict("axis", "x"));
	}

	int nx = this_img->get_xsize();

	if (maxshift < 0) {
		maxshift = nx / 10;
	}

	float angle_step =  params.set_default("angstep", 0.0f);
	if ( angle_step == 0 ) angle_step = atan2(2.0f, (float)nx);
	else {
		angle_step *= (float)EMConsts::deg2rad; //convert to radians
	}
	float trans_step =  params.set_default("transtep",1.0f);

	if (trans_step <= 0) throw InvalidParameterException("transstep must be greater than 0");
	if (angle_step <= 0) throw InvalidParameterException("angstep must be greater than 0");


	Dict shrinkfactor("n",2);
	EMData *this_img_shrink = this_img->process("math.medianshrink",shrinkfactor);
	EMData *to_shrunk = to->process("math.medianshrink",shrinkfactor);
	EMData *flipped_shrunk = flipped->process("math.medianshrink",shrinkfactor);

	int bestflip = 0;
	float bestdx = 0;
	float bestdy = 0;

	float bestang = 0;
	float bestval = FLT_MAX;

	int half_maxshift = maxshift / 2;


	for (int dy = -half_maxshift; dy <= half_maxshift; ++dy) {
		for (int dx = -half_maxshift; dx <= half_maxshift; ++dx) {
			if (hypot(dx, dy) <= maxshift) {
				for (float ang = -angle_step * 2.0f; ang <= (float)2 * M_PI; ang += angle_step * 4.0f) {
					EMData v(*this_img_shrink);
					Transform t(Dict("type","2d","alpha",static_cast<float>(ang*EMConsts::rad2deg)));
					t.set_trans((float)dx,(float)dy);
					v.transform(t);
// 					v.rotate_translate(ang*EMConsts::rad2deg, 0.0f, 0.0f, (float)dx, (float)dy, 0.0f);

					float lc = v.cmp(cmp_name, to_shrunk, cmp_params);

					if (lc < bestval) {
						bestval = lc;
						bestang = ang;
						bestdx = (float)dx;
						bestdy = (float)dy;
						bestflip = 0;
					}

					lc = v.cmp(cmp_name,flipped_shrunk , cmp_params);
					if (lc < bestval) {
						bestval = lc;
						bestang = ang;
						bestdx = (float)dx;
						bestdy = (float)dy;
						bestflip = 1;
					}
				}
			}
		}
	}

	if( to_shrunk )
	{
		delete to_shrunk;
		to_shrunk = 0;
	}
	if( flipped_shrunk )
	{
		delete flipped_shrunk;
		flipped_shrunk = 0;
	}
	if( this_img_shrink )
	{
		delete this_img_shrink;
		this_img_shrink = 0;
	}

	bestdx *= 2;
	bestdy *= 2;
	bestval = FLT_MAX;

	float bestdx2 = bestdx;
	float bestdy2 = bestdy;
	float bestang2 = bestang;

	for (float dy = bestdy2 - 3; dy <= bestdy2 + 3; dy += trans_step) {
		for (float dx = bestdx2 - 3; dx <= bestdx2 + 3; dx += trans_step) {
			if (hypot(dx, dy) <= maxshift) {
				for (float ang = bestang2 - angle_step * 6.0f; ang <= bestang2 + angle_step * 6.0f; ang += angle_step) {
					EMData v(*this_img);
					Transform t(Dict("type","2d","alpha",static_cast<float>(ang*EMConsts::rad2deg)));
					t.set_trans(dx,dy);
					v.transform(t);
// 					v.rotate_translate(ang*EMConsts::rad2deg, 0.0f, 0.0f, (float)dx, (float)dy, 0.0f);

					float lc = v.cmp(cmp_name, to, cmp_params);

					if (lc < bestval) {
						bestval = lc;
						bestang = ang;
						bestdx = dx;
						bestdy = dy;
						bestflip = 0;
					}

					lc = v.cmp(cmp_name, flipped, cmp_params);

					if (lc < bestval) {
						bestval = lc;
						bestang = ang;
						bestdx = dx;
						bestdy = dy;
						bestflip = 1;
					}
				}
			}
		}
	}

	if (delete_flipped) { delete flipped; flipped = 0; }

	bestang *= (float)EMConsts::rad2deg;
	Transform t(Dict("type","2d","alpha",(float)bestang));
	t.set_trans(bestdx,bestdy);

	if (bestflip) {
		t.set_mirror(true);
	}

	EMData* rslt = this_img->process("xform",Dict("transform",&t));
	rslt->set_attr("xform.align2d",&t);

	return rslt;
}

EMData* SymAlignProcessor::align(EMData * this_img, EMData *to, const string & cmp_name, const Dict& cmp_params) const
{

	// Set parms
	float dphi = params.set_default("dphi",10.f);
	float lphi = params.set_default("lphi",0.0f);
	float uphi = params.set_default("uphi",359.9f);

	Dict d;
	d["inc_mirror"] = true;
	d["delta"] = params.set_default("delta",10.f);

	//Genrate points on a sphere in an asymmetric unit
	Symmetry3D* sym = Factory<Symmetry3D>::get((string)params.set_default("sym","c1"));
	vector<Transform> transforms = sym->gen_orientations((string)params.set_default("orientgen","eman"),d);

	//Genrate symmetry related orritenations
	vector<Transform> syms = Symmetry3D::get_symmetries((string)params["sym"]);

	float bestquality = 0.0f;
	EMData* bestimage = 0;
	for(vector<Transform>::const_iterator trans_it = transforms.begin(); trans_it != transforms.end(); trans_it++) {
		Dict tparams = trans_it->get_params("eman");
		Transform t(tparams);
		for( float phi = lphi; phi < uphi; phi += dphi ) {
			tparams["phi"] = phi;
			t.set_rotation(tparams);

			//Get the averagaer
			Averager* imgavg = Factory<Averager>::get((string)params.set_default("avger","mean"));
			//Now make the averages
			for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {
				Transform sympos = (*it)*t;
				EMData* transformed = this_img->process("xform",Dict("transform",&sympos));
				imgavg->add_image(transformed);
				delete transformed;
			}

			EMData* symptcl=imgavg->finish();
			delete imgavg;
			//See which average is the best
			float quality = symptcl->get_attr("sigma");
			cout << quality << " " << phi << endl;
			if(quality > bestquality) {
				bestquality = quality;
				bestimage = symptcl;
			} else {
				delete symptcl;
			}
		}
	}
	if(sym != 0) delete sym;

	return bestimage;
}

static double refalifn(const gsl_vector * v, void *params)
{
	Dict *dict = (Dict *) params;

	double x = gsl_vector_get(v, 0);
	double y = gsl_vector_get(v, 1);
	double a = gsl_vector_get(v, 2);

	EMData *this_img = (*dict)["this"];
	EMData *with = (*dict)["with"];
	bool mirror = (*dict)["mirror"];

	Transform t(Dict("type","2d","alpha",static_cast<float>(a)));
	t.set_trans((float)x,(float)y);
	t.set_mirror(mirror);
	if (v->size>3) {
		float sca=(float)gsl_vector_get(v, 3);
		if (sca<.7 || sca>1.3) return 1.0e20;
		t.set_scale((float)gsl_vector_get(v, 3));
	}
	EMData *tmp = this_img->process("xform",Dict("transform",&t));
	if (dict->has_key("mask")) tmp->mult(*(EMData *)((*dict)["mask"]));

//	printf("GSL %f %f %f %d %f\n",x,y,a,mirror,(float)gsl_vector_get(v, 3));
	Cmp* c = (Cmp*) ((void*)(*dict)["cmp"]);
	double result = c->cmp(tmp,with);

	if (tmp != 0) delete tmp;

	return result;
}

static void refalidf(const gsl_vector * v, void *params,gsl_vector * df) {
	// we do this using a simple local difference estimate due to the expense of the calculation.
	// The step has to be large enough for the similarity metric
	// To provide an accurate change in value.
	static double lstep[4] = { 0.05, 0.05, 0.1, 0.01 };

	gsl_vector *vc = gsl_vector_alloc(v->size);
	gsl_vector_memcpy(vc,v);

	double f = refalifn(v,params);
	for (unsigned int i=0; i<v->size; i++) {
		// double *vp = gsl_vector_ptr(vc,i);
		double vp = gsl_vector_get(vc,i);
		// *vp+=lstep[i];
		gsl_vector_set(vc,i,vp+lstep[i]);
		double f2 = refalifn(vc,params);
		// *vp-=lstep[i];
		gsl_vector_set(vc,i,vp);

		gsl_vector_set(df,i,(f2-f)/lstep[i]);
	}

	gsl_vector_free(vc);
	return;
}

static void refalifdf(const gsl_vector * v, void *params, double * f, gsl_vector * df) {
	// we do this using a simple local difference estimate due to the expense of the calculation.
	// The step has to be large enough for the similarity metric
	// To provide an accurate change in value.
	static double lstep[4] = { 0.05, 0.05, 0.1, 0.01 };

	gsl_vector *vc = gsl_vector_alloc(v->size);
	gsl_vector_memcpy(vc,v);

	*f = refalifn(v,params);
	for (unsigned int i=0; i<v->size; i++) {
		// double *vp = gsl_vector_ptr(vc,i);
		double vp = gsl_vector_get(vc,i);
		// *vp+=lstep[i];
		gsl_vector_set(vc,i,vp+lstep[i]);
		double f2 = refalifn(vc,params);
		// *vp-=lstep[i];
		gsl_vector_set(vc,i,vp);

		gsl_vector_set(df,i,(f2-*f)/lstep[i]);
	}

	gsl_vector_free(vc);
	return;
}

static double refalifnfast(const gsl_vector * v, void *params)
{
	Dict *dict = (Dict *) params;
	EMData *this_img = (*dict)["this"];
	EMData *img_to = (*dict)["with"];
	bool mirror = (*dict)["mirror"];

	double x = gsl_vector_get(v, 0);
	double y = gsl_vector_get(v, 1);
	double a = gsl_vector_get(v, 2);

	double r = this_img->dot_rotate_translate(img_to, (float)x, (float)y, (float)a, mirror);
	int nsec = this_img->get_xsize() * this_img->get_ysize();
	double result = 1.0 - r / nsec;

// 	cout << result << " x " << x << " y " << y << " az " << a <<  endl;
	return result;
}

EMData *RefineAligner::align(EMData * this_img, EMData *to,
	const string & cmp_name, const Dict& cmp_params) const
{

	if (!to) {
		return 0;
	}

	EMData *result;
	int mode = params.set_default("mode", 0);
	float saz = 0.0;
	float sdx = 0.0;
	float sdy = 0.0;
	float sscale = 1.0;
	bool mirror = false;
	Transform* t;
	if (params.has_key("xform.align2d") ) {
		t = params["xform.align2d"];
		Dict params = t->get_params("2d");
		saz = params["alpha"];
		sdx = params["tx"];
		sdy = params["ty"];
		mirror = params["mirror"];
		sscale = params["scale"];
	} else {
		t = new Transform(); // is the identity
	}

	// We do this to prevent the GSL routine from crashing on an invalid alignment
	if ((float)(this_img->get_attr("sigma"))==0.0 || (float)(to->get_attr("sigma"))==0.0) {
		result = this_img->process("xform",Dict("transform",t));
		result->set_attr("xform.align2d",t);
		delete t;
		return result;
	}

	float stepx = params.set_default("stepx",1.0f);
	float stepy = params.set_default("stepy",1.0f);
	// Default step is 5 degree - note in EMAN1 it was 0.1 radians
	float stepaz = params.set_default("stepaz",5.0f);
	float stepscale = params.set_default("stepscale",0.0f);

	int np = 3;
	if (stepscale!=0.0) np++;
	Dict gsl_params;
	gsl_params["this"] = this_img;
	gsl_params["with"] = to;
	gsl_params["snr"]  = params["snr"];
	gsl_params["mirror"] = mirror;
	if (params.has_key("mask")) gsl_params["mask"]=params["mask"];

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
	gsl_vector *ss = gsl_vector_alloc(np);


	gsl_vector_set(ss, 0, stepx);
	gsl_vector_set(ss, 1, stepy);
	gsl_vector_set(ss, 2, stepaz);
	if (stepscale!=0.0) gsl_vector_set(ss,3,stepscale);

	gsl_vector *x = gsl_vector_alloc(np);
	gsl_vector_set(x, 0, sdx);
	gsl_vector_set(x, 1, sdy);
	gsl_vector_set(x, 2, saz);
	if (stepscale!=0.0) gsl_vector_set(x,3,1.0);

	Cmp *c = 0;

	gsl_multimin_function minex_func;
	if (mode == 2) {
		minex_func.f = &refalifnfast;
	}
	else {
		c = Factory < Cmp >::get(cmp_name, cmp_params);
		gsl_params["cmp"] = (void *) c;
		minex_func.f = &refalifn;
	}

	minex_func.n = np;
	minex_func.params = (void *) &gsl_params;

	gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, np);
	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	int rval = GSL_CONTINUE;
	int status = GSL_SUCCESS;
	int iter = 1;

	float precision = params.set_default("precision",0.04f);
	int maxiter = params.set_default("maxiter",28);

//	printf("Refine sx=%1.2f sy=%1.2f sa=%1.2f prec=%1.4f maxit=%d\n",stepx,stepy,stepaz,precision,maxiter);
//	printf("%1.2f %1.2f %1.1f  ->",(float)gsl_vector_get(s->x, 0),(float)gsl_vector_get(s->x, 1),(float)gsl_vector_get(s->x, 2));

	while (rval == GSL_CONTINUE && iter < maxiter) {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status) {
			break;
		}
		rval = gsl_multimin_test_size(gsl_multimin_fminimizer_size(s), precision);
	}

	int maxshift = params.set_default("maxshift",-1);

	if (maxshift <= 0) {
		maxshift = this_img->get_xsize() / 4;
	}
	float fmaxshift = static_cast<float>(maxshift);
	if ( fmaxshift >= fabs((float)gsl_vector_get(s->x, 0)) && fmaxshift >= fabs((float)gsl_vector_get(s->x, 1)) && (stepscale==0 || (((float)gsl_vector_get(s->x, 3))<1.3 && ((float)gsl_vector_get(s->x, 3))<0.7))  )
	{
//		printf(" Refine good %1.2f %1.2f %1.1f\n",(float)gsl_vector_get(s->x, 0),(float)gsl_vector_get(s->x, 1),(float)gsl_vector_get(s->x, 2));
		Transform  tsoln(Dict("type","2d","alpha",(float)gsl_vector_get(s->x, 2)));
		tsoln.set_mirror(mirror);
		tsoln.set_trans((float)gsl_vector_get(s->x, 0),(float)gsl_vector_get(s->x, 1));
		if (stepscale!=0.0) tsoln.set_scale((float)gsl_vector_get(s->x, 3));
		result = this_img->process("xform",Dict("transform",&tsoln));
		result->set_attr("xform.align2d",&tsoln);
	} else { // The refine aligner failed - this shift went beyond the max shift
//		printf(" Refine Failed %1.2f %1.2f %1.1f\n",(float)gsl_vector_get(s->x, 0),(float)gsl_vector_get(s->x, 1),(float)gsl_vector_get(s->x, 2));
		result = this_img->process("xform",Dict("transform",t));
		result->set_attr("xform.align2d",t);
	}

	delete t;
	t = 0;

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);

	if (c != 0) delete c;
	return result;
}

EMData *RefineAlignerCG::align(EMData * this_img, EMData *to,
	const string & cmp_name, const Dict& cmp_params) const
{

	if (!to) {
		return 0;
	}

	EMData *result;
	int mode = params.set_default("mode", 0);
	float saz = 0.0;
	float sdx = 0.0;
	float sdy = 0.0;
	float sscale = 1.0;
	bool mirror = false;
	Transform* t;
	if (params.has_key("xform.align2d") ) {
		t = params["xform.align2d"];
		Dict params = t->get_params("2d");
		saz = params["alpha"];
		sdx = params["tx"];
		sdy = params["ty"];
		mirror = params["mirror"];
		sscale = params["scale"];
	} else {
		t = new Transform(); // is the identity
	}

	// We do this to prevent the GSL routine from crashing on an invalid alignment
	if ((float)(this_img->get_attr("sigma"))==0.0 || (float)(to->get_attr("sigma"))==0.0) {
		result = this_img->process("xform",Dict("transform",t));
		result->set_attr("xform.align2d",t);
		delete t;
		return result;
	}

	float step = params.set_default("step",0.1f);
	float stepscale = params.set_default("stepscale",0.0f);

	int np = 3;
	if (stepscale!=0.0) np++;
	Dict gsl_params;
	gsl_params["this"] = this_img;
	gsl_params["with"] = to;
	gsl_params["snr"]  = params["snr"];
	gsl_params["mirror"] = mirror;
	if (params.has_key("mask")) gsl_params["mask"]=params["mask"];

	const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs;

	gsl_vector *x = gsl_vector_alloc(np);
	gsl_vector_set(x, 0, sdx);
	gsl_vector_set(x, 1, sdy);
	gsl_vector_set(x, 2, saz);
	if (stepscale!=0.0) gsl_vector_set(x,3,1.0);

	Cmp *c = 0;

	gsl_multimin_function_fdf minex_func;
	if (mode == 2) {
		minex_func.f = &refalifnfast;
	}
	else {
		c = Factory < Cmp >::get(cmp_name, cmp_params);
		gsl_params["cmp"] = (void *) c;
		minex_func.f = &refalifn;
	}

	minex_func.df = &refalidf;
	minex_func.fdf = &refalifdf;
	minex_func.n = np;
	minex_func.params = (void *) &gsl_params;

	gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, np);
	gsl_multimin_fdfminimizer_set(s, &minex_func, x, step, 0.001f);

	int rval = GSL_CONTINUE;
	int status = GSL_SUCCESS;
	int iter = 1;

	float precision = params.set_default("precision",0.02f);
	int maxiter = params.set_default("maxiter",12);
	int verbose = params.set_default("verbose",0);

//	printf("Refine sx=%1.2f sy=%1.2f sa=%1.2f prec=%1.4f maxit=%d\n",stepx,stepy,stepaz,precision,maxiter);
//	printf("%1.2f %1.2f %1.1f  ->",(float)gsl_vector_get(s->x, 0),(float)gsl_vector_get(s->x, 1),(float)gsl_vector_get(s->x, 2));

	while (rval == GSL_CONTINUE && iter < maxiter) {
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);
		if (status) {
			break;
		}
		rval = gsl_multimin_test_gradient (s->gradient, precision);
//		if (verbose>2) printf("GSL %d. %1.3f %1.3f %1.3f   %1.3f\n",iter,gsl_vector_get(s->x,0),gsl_vector_get(s->x,1),gsl_vector_get(s->x,2),s->gradient[0]);
	}

	int maxshift = params.set_default("maxshift",-1);

	if (maxshift <= 0) {
		maxshift = this_img->get_xsize() / 4;
	}
	float fmaxshift = static_cast<float>(maxshift);
	if ( fmaxshift >= fabs((float)gsl_vector_get(s->x, 0)) && fmaxshift >= fabs((float)gsl_vector_get(s->x, 1)) && (stepscale==0 || (((float)gsl_vector_get(s->x, 3))<1.3 && ((float)gsl_vector_get(s->x, 3))<0.7))  )
	{
		if (verbose>0) printf(" Refine good (%d) %1.2f %1.2f %1.1f\n",iter,(float)gsl_vector_get(s->x, 0),(float)gsl_vector_get(s->x, 1),(float)gsl_vector_get(s->x, 2));
		Transform  tsoln(Dict("type","2d","alpha",(float)gsl_vector_get(s->x, 2)));
		tsoln.set_mirror(mirror);
		tsoln.set_trans((float)gsl_vector_get(s->x, 0),(float)gsl_vector_get(s->x, 1));
		if (stepscale!=0.0) tsoln.set_scale((float)gsl_vector_get(s->x, 3));
		result = this_img->process("xform",Dict("transform",&tsoln));
		result->set_attr("xform.align2d",&tsoln);
	} else { // The refine aligner failed - this shift went beyond the max shift
		if (verbose>1) printf(" Refine Failed %1.2f %1.2f %1.1f\n",(float)gsl_vector_get(s->x, 0),(float)gsl_vector_get(s->x, 1),(float)gsl_vector_get(s->x, 2));
		result = this_img->process("xform",Dict("transform",t));
		result->set_attr("xform.align2d",t);
	}

	delete t;
	t = 0;

	gsl_vector_free(x);
	gsl_multimin_fdfminimizer_free(s);

	if (c != 0) delete c;
	return result;
}

static Transform refalin3d_perturbquat(const Transform*const t, const float& spincoeff, const float& n0, const float& n1, const float& n2, const float& x, const float& y, const float& z)
{
	Vec3f normal(n0,n1,n2);
	normal.normalize();

	float omega = spincoeff*sqrt(n0*n0 + n1*n1 + n2*n2); // Here we compute the spin by the rotation axis vector length
	Dict d;
	d["type"] = "spin";
	d["omega"] = omega;
	d["n1"] = normal[0];
	d["n2"] = normal[1];
	d["n3"] = normal[2];
	//cout << omega << " " << normal[0] << " " << normal[1] << " " << normal[2] << " " << n0 << " " << n1 << " " << n2 << endl;

	Transform q(d);
	q.set_trans((float)x,(float)y,(float)z);

	q = q*(*t); //compose transforms

	return q;
}

static double symquat(const gsl_vector * v, void *params)
{
	Dict *dict = (Dict *) params;

 	double n0 = gsl_vector_get(v, 0);
 	double n1 = gsl_vector_get(v, 1);
	double n2 = gsl_vector_get(v, 2);
	double x = gsl_vector_get(v, 3);
	double y = gsl_vector_get(v, 4);
	double z = gsl_vector_get(v, 5);

	EMData* volume = (*dict)["volume"];
	float spincoeff = (*dict)["spincoeff"];
	Transform* t = (*dict)["transform"];

	Transform soln = refalin3d_perturbquat(t,spincoeff,(float)n0,(float)n1,(float)n2,(float)x,(float)y,(float)z);

	EMData *tmp = volume->process("xform",Dict("transform",&soln));
	EMData *symtmp = tmp->process("xform.applysym",Dict("sym",(*dict)["sym"]));
	Cmp* c = (Cmp*) ((void*)(*dict)["cmp"]);
	double result = c->cmp(symtmp,tmp);
	delete tmp;
	delete symtmp;

	//cout << result << endl;
	return result;
}

static double refalifn3dquat(const gsl_vector * v, void *params)
{
	Dict *dict = (Dict *) params;

 	double n0 = gsl_vector_get(v, 0);
 	double n1 = gsl_vector_get(v, 1);
	double n2 = gsl_vector_get(v, 2);
	double x = gsl_vector_get(v, 3);
	double y = gsl_vector_get(v, 4);
	double z = gsl_vector_get(v, 5);

	EMData *this_img = (*dict)["this"];
	EMData *with = (*dict)["with"];

	Transform* t = (*dict)["transform"];
	float spincoeff = (*dict)["spincoeff"];

	Transform soln = refalin3d_perturbquat(t,spincoeff,(float)n0,(float)n1,(float)n2,(float)x,(float)y,(float)z);

	EMData *tmp = this_img->process("xform",Dict("transform",&soln));
	Cmp* c = (Cmp*) ((void*)(*dict)["cmp"]);
	double result = c->cmp(tmp,with);
	if ( tmp != 0 ) delete tmp;

	//cout << result << endl;
	return result;
}

EMData* SymAlignProcessorQuat::align(EMData * volume, EMData *to, const string & cmp_name, const Dict& cmp_params) const
{
	//Get pretransform
	Transform* t;
	if (params.has_key("xform.align3d") ) {
		t = params["xform.align3d"];
	}else {
		t = new Transform(); // is the identity
	}

	float sdi = 0.0;
	float sdj = 0.0;
	float sdk = 0.0;
	float sdx = 0.0;
	float sdy = 0.0;
	float sdz = 0.0;

	float spincoeff =  params.set_default("spin_coeff",10.0f); // spin coefficient, controls speed of convergence (sort of)

	int np = 6; // the number of dimensions
	Dict gsl_params;
	gsl_params["volume"] = volume;
	gsl_params["transform"] = t;
	gsl_params["sym"] = params.set_default("sym","c1");
	gsl_params["spincoeff"] = spincoeff;

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
	gsl_vector *ss = gsl_vector_alloc(np);

	float stepi = params.set_default("stepn0",1.0f); // doesn't really matter b/c the vecor part will be normalized anyway
	float stepj = params.set_default("stepn1",1.0f); // doesn't really matter b/c the vecor part will be normalized anyway
	float stepk = params.set_default("stepn2",1.0f); // doesn't really matter b/c the vecor part will be normalized anyway
	float stepx = params.set_default("stepx",1.0f);
	float stepy = params.set_default("stepy",1.0f);
	float stepz = params.set_default("stepz",1.0f);

	gsl_vector_set(ss, 0, stepi);
	gsl_vector_set(ss, 1, stepj);
	gsl_vector_set(ss, 2, stepk);
	gsl_vector_set(ss, 3, stepx);
	gsl_vector_set(ss, 4, stepy);
	gsl_vector_set(ss, 5, stepz);

	gsl_vector *x = gsl_vector_alloc(np);
	gsl_vector_set(x, 0, sdi);
	gsl_vector_set(x, 1, sdj);
	gsl_vector_set(x, 2, sdk);
	gsl_vector_set(x, 3, sdx);
	gsl_vector_set(x, 4, sdy);
	gsl_vector_set(x, 5, sdz);

	gsl_multimin_function minex_func;
	Cmp *c = Factory < Cmp >::get(cmp_name, cmp_params);
	gsl_params["cmp"] = (void *) c;
	minex_func.f = &symquat;
	minex_func.n = np;
	minex_func.params = (void *) &gsl_params;
	gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, np);
	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	int rval = GSL_CONTINUE;
	int status = GSL_SUCCESS;
	int iter = 1;

	float precision = params.set_default("precision",0.01f);
	int maxiter = params.set_default("maxiter",100);
	while (rval == GSL_CONTINUE && iter < maxiter) {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status) {
			break;
		}
		rval = gsl_multimin_test_size(gsl_multimin_fminimizer_size(s), precision);
	}

	int maxshift = params.set_default("maxshift",-1);

	if (maxshift <= 0) {
		maxshift = volume->get_xsize() / 4;
	}
	float fmaxshift = static_cast<float>(maxshift);

	EMData *result;
	if ( fmaxshift >= (float)gsl_vector_get(s->x, 0) && fmaxshift >= (float)gsl_vector_get(s->x, 1)  && fmaxshift >= (float)gsl_vector_get(s->x, 2))
	{
		float n0 = (float)gsl_vector_get(s->x, 0);
		float n1 = (float)gsl_vector_get(s->x, 1);
		float n2 = (float)gsl_vector_get(s->x, 2);
		float x = (float)gsl_vector_get(s->x, 3);
		float y = (float)gsl_vector_get(s->x, 4);
		float z = (float)gsl_vector_get(s->x, 5);

		Transform tsoln = refalin3d_perturbquat(t,spincoeff,n0,n1,n2,x,y,z);

		result = volume->process("xform",Dict("transform",&tsoln));
		result->set_attr("xform.align3d",&tsoln);
		EMData *tmpsym = result->process("xform.applysym",Dict("sym",gsl_params["sym"]));
		result->set_attr("score", result->cmp(cmp_name,tmpsym,cmp_params));
		delete tmpsym;
	} else { // The refine aligner failed - this shift went beyond the max shift
		result = volume->process("xform",Dict("transform",t));
		result->set_attr("xform.align3d",t);
		result->set_attr("score",0.0);
	}

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);

	if (c != 0) delete c;
	delete t;

	return result;
}

EMData* Refine3DAlignerQuaternion::align(EMData * this_img, EMData *to,
	const string & cmp_name, const Dict& cmp_params) const
{

	if (!to || !this_img) throw NullPointerException("Input image is null"); // not sure if this is necessary, it was there before I started

	if (to->get_ndim() != 3 || this_img->get_ndim() != 3) throw ImageDimensionException("The Refine3D aligner only works for 3D images");

#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		if(!this_img->getcudarwdata()) this_img->copy_to_cuda();
		if(!to->getcudarwdata()) to->copy_to_cuda();
	}
#endif

	float sdi = 0.0;
	float sdj = 0.0;
	float sdk = 0.0;
	float sdx = 0.0;
	float sdy = 0.0;
	float sdz = 0.0;
	bool mirror = false;

	Transform* t;
	if (params.has_key("xform.align3d") ) {
		// Unlike the 2d refine aligner, this class doesn't require the starting transform's
		// parameters to form the starting guess. Instead the Transform itself
		// is perturbed carefully (using quaternion rotation) to overcome problems that arise
		// when you use orthogonally-based Euler angles
		t = params["xform.align3d"];
	}else {
		t = new Transform(); // is the identity
	}

	float spincoeff =  params.set_default("spin_coeff",10.0f); // spin coefficient, controls speed of convergence (sort of)

	int np = 6; // the number of dimensions
	Dict gsl_params;
	gsl_params["this"] = this_img;
	gsl_params["with"] = to;
	gsl_params["snr"]  = params["snr"];
	gsl_params["mirror"] = mirror;
	gsl_params["transform"] = t;
	gsl_params["spincoeff"] = spincoeff;
	Dict altered_cmp_params(cmp_params);

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
	gsl_vector *ss = gsl_vector_alloc(np);

	float stepi = params.set_default("stepn0",1.0f); // doesn't really matter b/c the vecor part will be normalized anyway
	float stepj = params.set_default("stepn1",1.0f); // doesn't really matter b/c the vecor part will be normalized anyway
	float stepk = params.set_default("stepn2",1.0f); // doesn't really matter b/c the vecor part will be normalized anyway
	float stepx = params.set_default("stepx",1.0f);
	float stepy = params.set_default("stepy",1.0f);
	float stepz = params.set_default("stepz",1.0f);

	//gsl_vector_set(ss, 0, stepw);
	gsl_vector_set(ss, 0, stepi);
	gsl_vector_set(ss, 1, stepj);
	gsl_vector_set(ss, 2, stepk);
	gsl_vector_set(ss, 3, stepx);
	gsl_vector_set(ss, 4, stepy);
	gsl_vector_set(ss, 5, stepz);

	gsl_vector *x = gsl_vector_alloc(np);
	gsl_vector_set(x, 0, sdi);
	gsl_vector_set(x, 1, sdj);
	gsl_vector_set(x, 2, sdk);
	gsl_vector_set(x, 3, sdx);
	gsl_vector_set(x, 4, sdy);
	gsl_vector_set(x, 5, sdz);

	gsl_multimin_function minex_func;
	Cmp *c = Factory < Cmp >::get(cmp_name, altered_cmp_params);

	gsl_params["cmp"] = (void *) c;
	minex_func.f = &refalifn3dquat;

	minex_func.n = np;
	minex_func.params = (void *) &gsl_params;

	gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, np);
	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	int rval = GSL_CONTINUE;
	int status = GSL_SUCCESS;
	int iter = 1;

	float precision = params.set_default("precision",0.01f);
	int maxiter = params.set_default("maxiter",100);
	while (rval == GSL_CONTINUE && iter < maxiter) {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status) {
			break;
		}
		rval = gsl_multimin_test_size(gsl_multimin_fminimizer_size(s), precision);
	}

	int maxshift = params.set_default("maxshift",-1);

	if (maxshift <= 0) {
		maxshift = this_img->get_xsize() / 4;
	}
	float fmaxshift = static_cast<float>(maxshift);

	EMData *result;
	if ( fmaxshift >= (float)gsl_vector_get(s->x, 0) && fmaxshift >= (float)gsl_vector_get(s->x, 1)  && fmaxshift >= (float)gsl_vector_get(s->x, 2))
	{
		float n0 = (float)gsl_vector_get(s->x, 0);
		float n1 = (float)gsl_vector_get(s->x, 1);
		float n2 = (float)gsl_vector_get(s->x, 2);
		float x = (float)gsl_vector_get(s->x, 3);
		float y = (float)gsl_vector_get(s->x, 4);
		float z = (float)gsl_vector_get(s->x, 5);

		Transform tsoln = refalin3d_perturbquat(t,spincoeff,n0,n1,n2,x,y,z);

		result = this_img->process("xform",Dict("transform",&tsoln));
		result->set_attr("xform.align3d",&tsoln);
		result->set_attr("score", result->cmp(cmp_name,to,cmp_params));

	 //coda goes here
	} else { // The refine aligner failed - this shift went beyond the max shift
		result = this_img->process("xform",Dict("transform",t));
		result->set_attr("xform.align3d",t);
		result->set_attr("score",0.0);
	}

	//EMData *result = this_img->process("xform",Dict("transform",t));
	delete t;
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);

	if (c != 0) delete c;

	return result;
}

EMData* Refine3DAlignerGrid::align(EMData * this_img, EMData *to,
	const string & cmp_name, const Dict& cmp_params) const
{
	if ( this_img->get_ndim() != 3 || to->get_ndim() != 3 ) {
		throw ImageDimensionException("This aligner only works for 3D images");
	}

	Transform* t;
	if (params.has_key("xform.align3d") ) {
		// Unlike the 2d refine aligner, this class doesn't require the starting transform's
		// parameters to form the starting guess. Instead the Transform itself
		// is perturbed carefully (using quaternion rotation) to overcome problems that arise
		// when you use orthogonally-based Euler angles
		t = params["xform.align3d"];
	}else {
		t = new Transform(); // is the identity
	}

	int searchx = 0;
	int searchy = 0;
	int searchz = 0;
	bool dotrans = params.set_default("dotrans",1);
	if (params.has_key("search")) {
		vector<string> check;
		check.push_back("searchx");
		check.push_back("searchy");
		check.push_back("searchz");
		for(vector<string>::const_iterator cit = check.begin(); cit != check.end(); ++cit) {
			if (params.has_key(*cit)) throw InvalidParameterException("The search parameter is mutually exclusive of the searchx, searchy, and searchz parameters");
		}
		int search  = params["search"];
		searchx = search;
		searchy = search;
		searchz = search;
	} else {
		searchx = params.set_default("searchx",3);
		searchy = params.set_default("searchy",3);
		searchz = params.set_default("searchz",3);
	}

	float delta = params.set_default("delta",1.0f);
	float range = params.set_default("range",10.0f);
	bool verbose = params.set_default("verbose",false);

	bool tomography = (cmp_name == "ccc.tomo") ? 1 : 0;
	EMData * tofft = 0;
	if(dotrans || tomography){
		tofft = to->do_fft();
	}

#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		if(!this_img->getcudarodata()) this_img->copy_to_cudaro(); // This is safer
		if(!to->getcudarwdata()) to->copy_to_cuda();
		if(to->getcudarwdata()){if(tofft) tofft->copy_to_cuda();}
	}
#endif

	Dict d;
	d["type"] = "eman"; // d is used in the loop below
	Dict best;
//	best["score"] = numeric_limits<float>::infinity();
	best["score"] = 1.0e37;
	bool use_cpu = true;
	Transform tran = Transform();
	Cmp* c = Factory <Cmp>::get(cmp_name, cmp_params);

	for ( float alt = 0; alt < range; alt += delta) {
		// now compute a sane az step size
		float saz = 360;
		if(alt != 0) saz = delta/sin(alt*M_PI/180.0f); // This gives consistent az step sizes(arc lengths)
		for ( float az = 0; az < 360; az += saz ){
			if (verbose) {
				cout << "Trying angle alt " << alt << " az " << az << endl;
			}
			// account for any changes in az
			for( float phi = -range-az; phi < range-az; phi += delta ) {
				d["alt"] = alt;
				d["phi"] = phi;
				d["az"] = az;
				Transform tr(d);
				tr = tr*(*t);	// compose transforms, this moves to the pole (aprox)

				//EMData* transformed = this_img->process("xform",Dict("transform",&tr));
				EMData* transformed;
				transformed = this_img->process("xform",Dict("transform",&tr));

				//need to do things a bit diffrent if we want to compare two tomos
				float score = 0.0f;
				if(dotrans || tomography){
					EMData* ccf = transformed->calc_ccf(tofft);
#ifdef EMAN2_USING_CUDA
					if(EMData::usecuda == 1){
						use_cpu = false;
						CudaPeakInfo* data = calc_max_location_wrap_cuda(ccf->getcudarwdata(), ccf->get_xsize(), ccf->get_ysize(), ccf->get_zsize(), searchx, searchy, searchz);
						tran.set_trans((float)-data->px, (float)-data->py, (float)-data->pz);
						//CudaPeakInfoFloat* data = calc_max_location_wrap_intp_cuda(ccf->getcudarwdata(), ccf->get_xsize(), ccf->get_ysize(), ccf->get_zsize(), searchx, searchy, searchz);
						//tran.set_trans(-data->xintp, -data->yintp, -data->zintp);
						tr = tran*tr; // to reflect the fact that we have done a rotation first and THEN a transformation
						if (tomography) {
							float2 stats = get_stats_cuda(ccf->getcudarwdata(), ccf->get_xsize(), ccf->get_ysize(), ccf->get_zsize());
							score = -(data->peak - stats.x)/sqrt(stats.y); // Normalize, this is better than calling the norm processor since we only need to normalize one point
						} else {
							score = -data->peak;
						}
						delete data;
					}
#endif
					if(use_cpu){
						if(tomography) ccf->process_inplace("normalize");
						//vector<float> fpoint = ccf->calc_max_location_wrap_intp(searchx,searchy,searchz);
						//tran.set_trans(-fpoint[0], -fpoint[1], -fpoint[2]);
						//score = -fpoint[3];
						IntPoint point = ccf->calc_max_location_wrap(searchx,searchy,searchz);
						tran.set_trans((float)-point[0], (float)-point[1], (float)-point[2]);
						score = -ccf->get_value_at_wrap(point[0], point[1], point[2]);
						tr = tran*tr;// to reflect the fact that we have done a rotation first and THEN a transformation

					}
					delete ccf; ccf =0;
					delete transformed; transformed = 0;// this is to stop a mem leak
				}

				if(!tomography){
					if(!transformed) transformed = this_img->process("xform",Dict("transform",&tr)); // we are returning a map
					score = c->cmp(to,transformed);
					delete transformed; transformed = 0;// this is to stop a mem leak
				}

				if(score < float(best["score"])) {
//					printf("%f\n",score);
					best["score"] = score;
					best["xform.align3d"] = &tr; // I wonder if this will cause a mem leak?
				}
			}
		}
	}

	if(tofft) {delete tofft; tofft = 0;}
	if (c != 0) delete c;

	//make aligned map;
	EMData* best_match = this_img->process("xform",Dict("transform", best["xform.align3d"])); // we are returning a map
	best_match->set_attr("xform.align3d", best["xform.align3d"]);
	best_match->set_attr("score", float(best["score"]));

	return best_match;

}

EMData* RT3DGridAligner::align(EMData * this_img, EMData *to, const string & cmp_name, const Dict& cmp_params) const
{

	vector<Dict> alis = xform_align_nbest(this_img,to,1,cmp_name,cmp_params);

	Dict t;
	Transform* tr = (Transform*) alis[0]["xform.align3d"];
	t["transform"] = tr;
	EMData* soln = this_img->process("xform",t);
	soln->set_attr("xform.align3d",tr);
	delete tr; tr = 0;

	return soln;

}

vector<Dict> RT3DGridAligner::xform_align_nbest(EMData * this_img, EMData * to, const unsigned int nsoln, const string & cmp_name, const Dict& cmp_params) const {

	if ( this_img->get_ndim() != 3 || to->get_ndim() != 3 ) {
		throw ImageDimensionException("This aligner only works for 3D images");
	}

	int searchx = 0;
	int searchy = 0;
	int searchz = 0;

	bool dotrans = params.set_default("dotrans",1);
	if (params.has_key("search")) {
		vector<string> check;
		check.push_back("searchx");
		check.push_back("searchy");
		check.push_back("searchz");
		for(vector<string>::const_iterator cit = check.begin(); cit != check.end(); ++cit) {
			if (params.has_key(*cit)) throw InvalidParameterException("The search parameter is mutually exclusive of the searchx, searchy, and searchz parameters");
		}
		int search  = params["search"];
		searchx = search;
		searchy = search;
		searchz = search;
	} else {
		searchx = params.set_default("searchx",3);
		searchy = params.set_default("searchy",3);
		searchz = params.set_default("searchz",3);
	}

	Transform* initxform;
	if (params.has_key("initxform") ) {
		// Unlike the 2d refine aligner, this class doesn't require the starting transform's
		// parameters to form the starting guess. Instead the Transform itself
		// is perturbed carefully (using quaternion rotation) to overcome problems that arise
		// when you use orthogonally-based Euler angles
		initxform = params["initxform"];
	}else {
		initxform = new Transform(); // is the identity
	}

	float lalt = params.set_default("alt0",0.0f);
	float laz = params.set_default("az0",0.0f);
	float lphi = params.set_default("phi0",0.0f);
	float ualt = params.set_default("alt1",180.0f); // I am using 179.9 rather than 180 to avoid resampling
	float uphi = params.set_default("phi1",360.0f); // I am using 359.9 rather than 180 to avoid resampling 0 = 360 (for perodic functions)
	float uaz = params.set_default("az1",360.0f);   // I am using 359.9 rather than 180 to avoid resampling 0 = 360 (for perodic functions)
	float dalt = params.set_default("dalt",10.f);
	float daz = params.set_default("daz",10.f);
	float dphi = params.set_default("dphi",10.f);
	bool verbose = params.set_default("verbose",false);

	//in case we arre aligning tomos
	Dict altered_cmp_params(cmp_params);
	if (cmp_name == "ccc.tomo") {
                altered_cmp_params.set_default("searchx", searchx);
		altered_cmp_params.set_default("searchy", searchy);
		altered_cmp_params.set_default("searchz", searchz);
		altered_cmp_params.set_default("norm", true);
	}

	vector<Dict> solns;
	if (nsoln == 0) return solns; // What was the user thinking?
	for (unsigned int i = 0; i < nsoln; ++i ) {
		Dict d;
		d["score"] = 1.e24;
		Transform t; // identity by default
		d["xform.align3d"] = &t; // deep copy is going on here
		solns.push_back(d);
	}

	bool tomography = (cmp_name == "ccc.tomo") ? 1 : 0;
	EMData * tofft = 0;
	if(dotrans || tomography){
		tofft = to->do_fft();
	}

#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		if(!this_img->getcudarodata()) this_img->copy_to_cudaro();  // safer call
		if(!to->getcudarwdata()) to->copy_to_cuda();
		if(to->getcudarwdata()){if(tofft) tofft->copy_to_cuda();}
	}
#endif

	Dict d;
	d["type"] = "eman"; // d is used in the loop below
	Transform trans = Transform();
	Cmp* c = Factory <Cmp>::get(cmp_name, cmp_params);
	bool use_cpu = true;
	for ( float alt = lalt; alt <= ualt; alt += dalt) {
		// An optimization for the range of az is made at the top of the sphere
		// If you think about it, this is just a coarse way of making this approach slightly more efficient
		for ( float az = laz; az < uaz; az += daz ){
			if (verbose) {
				cout << "Trying angle alt " << alt << " az " << az << endl;
			}
			for( float phi = lphi; phi < uphi; phi += dphi ) {
				d["alt"] = alt;
				d["phi"] = phi;
				d["az"] = az;
				Transform t(d);
				t = t*(*initxform);
				EMData* transformed = this_img->process("xform",Dict("transform",&t));

				//need to do things a bit diffrent if we want to compare two tomos
				float best_score = 0.0f;
				if(dotrans || tomography){
					EMData* ccf = transformed->calc_ccf(tofft);
#ifdef EMAN2_USING_CUDA
					if(EMData::usecuda == 1){
						use_cpu = false;;
						CudaPeakInfo* data = calc_max_location_wrap_cuda(ccf->getcudarwdata(), ccf->get_xsize(), ccf->get_ysize(), ccf->get_zsize(), searchx, searchy, searchz);
						trans.set_trans((float)-data->px, (float)-data->py, (float)-data->pz);
						t = trans*t;	//composite transfrom to reflect the fact that we have done a rotation first and THEN a transformation
						if (tomography) {
							float2 stats = get_stats_cuda(ccf->getcudarwdata(), ccf->get_xsize(), ccf->get_ysize(), ccf->get_zsize());
							best_score = -(data->peak - stats.x)/sqrt(stats.y); // Normalize, this is better than calling the norm processor since we only need to normalize one point
						} else {
							best_score = -data->peak;
						}
						delete data;
					}
#endif
					if(use_cpu){
						if(tomography) ccf->process_inplace("normalize");
						IntPoint point = ccf->calc_max_location_wrap(searchx,searchy,searchz);
						trans.set_trans((float)-point[0], (float)-point[1], (float)-point[2]);
						t = trans*t;	//composite transfrom to reflect the fact that we have done a rotation first and THEN a transformation
						best_score = -ccf->get_value_at_wrap(point[0], point[1], point[2]);
					}
					delete ccf; ccf =0;
					delete transformed; transformed = 0;
				}

				if(!tomography){
					if(!transformed) transformed = this_img->process("xform",Dict("transform",&t));
					best_score = c->cmp(to,transformed);
					delete transformed; transformed = 0;
				}

				unsigned int j = 0;
				for ( vector<Dict>::iterator it = solns.begin(); it != solns.end(); ++it, ++j ) {
					if ( (float)(*it)["score"] > best_score ) {  // Note greater than - EMAN2 preferes minimums as a matter of policy
						vector<Dict>::reverse_iterator rit = solns.rbegin();
						copy(rit+1,solns.rend()-j,rit);
						Dict& d = (*it);
						d["score"] = best_score;
						d["xform.align3d"] = &t;
						break;
					}
				}
			}
		}
	}

	if(tofft) {delete tofft; tofft = 0;}
	if (c != 0) delete c;

	return solns;

}
EMData* RT3DTreeAligner::align(EMData * this_img, EMData *to, const string & cmp_name, const Dict& cmp_params) const
{

 	vector<Dict> alis = xform_align_nbest(this_img,to,8,cmp_name,cmp_params);

 	Dict t;
 	Transform* tr = (Transform*) alis[0]["xform.align3d"];
 	t["transform"] = tr;
 	EMData* soln = this_img->process("xform",t);
 	soln->set_attr("xform.align3d",tr);

	return soln;

}

// NOTE - if symmetry is applied, it is critical that "to" be the volume which is already aligned to the symmetry axes (ie - the reference)
vector<Dict> RT3DTreeAligner::xform_align_nbest(EMData * this_img, EMData * to, const unsigned int nrsoln, const string & cmp_name, const Dict& cmp_params) const {
	if (nrsoln == 0) throw InvalidParameterException("ERROR (RT3DTreeAligner): nsoln must be >0"); // What was the user thinking?

	int nsoln = nrsoln*2;
	if (nrsoln<16) nsoln=32;		// we start with at least 32 solutions, but then gradually decrease with increasing scale
	
	// !!!!!! IMPORTANT NOTE - we are inverting the order of this and to here to match convention in other aligners, to compensate
	// the Transform is inverted before being returned
	EMData *base_this;
	EMData *base_to;
	if (this_img->is_complex()) base_to=this_img->copy();
	else {
		base_to=this_img->do_fft();
		base_to->process_inplace("xform.phaseorigin.tocorner");
	}
	
	if (to->is_complex()) base_this=to->copy();
	else {
		base_this=to->do_fft();
		base_this->process_inplace("xform.phaseorigin.tocorner");
	}
	
	float sigmathis = params.set_default("sigmathis",0.01f);
	float sigmato = params.set_default("sigmato",0.01f);
	int verbose = params.set_default("verbose",0);


	if (base_this->get_xsize()!=base_this->get_ysize()+2 || base_this->get_ysize()!=base_this->get_zsize()
		|| base_to->get_xsize()!=base_to->get_ysize()+2 || base_to->get_ysize()!=base_to->get_zsize()) throw InvalidCallException("ERROR (RT3DTreeAligner): requires cubic images with even numbered box sizes");

	base_this->process_inplace("xform.fourierorigin.tocenter");		// easier to chop out Fourier subvolumes
	base_to->process_inplace("xform.fourierorigin.tocenter");

	float apix=(float)this_img->get_attr("apix_x");
	int ny=this_img->get_ysize();

//	int downsample=floor(ny/20);		// Minimum shrunken box size is 20^3

	vector<float> s_score(nsoln,0.0f);
	vector<float> s_coverage(nsoln,0.0f);
	vector<float> s_step(nsoln*3,7.5f);
	vector<Transform> s_xform(nsoln);
	if (verbose>0) printf("%d solutions\n",nsoln);

	
//	float dstep[3] = {7.5,7.5,7.5};		// we take  steps for each of the 3 angles, may be positive or negative
	char *axname[] = {"az","alt","phi"};

	// We start with 32^3, 64^3 ...
	for (int sexp=4; sexp<10; sexp++) {
		int ss=pow(2.0,sexp);
		if (ss==16) ss=24;		// 16 may be too small, but 32 takes too long...
		if (ss>ny) ss=ny;
		if (verbose>0) printf("\nSize %d\n",ss);

		//ss=good_size(ny/ds);
		EMData *small_this=base_this->get_clip(Region(0,(ny-ss)/2,(ny-ss)/2,ss+2,ss,ss));
		EMData *small_to=  base_to->  get_clip(Region(0,(ny-ss)/2,(ny-ss)/2,ss+2,ss,ss));
		small_this->process_inplace("xform.fourierorigin.tocorner");					// after clipping back to canonical form
		small_this->process_inplace("filter.highpass.gauss",Dict("cutoff_pixels",4));
		small_this->process_inplace("filter.lowpass.gauss",Dict("cutoff_abs",0.33f));
		small_to->process_inplace("xform.fourierorigin.tocorner");
		small_to->process_inplace("filter.highpass.gauss",Dict("cutoff_pixels",4));
		small_to->process_inplace("filter.lowpass.gauss",Dict("cutoff_abs",0.33f));

		// these are cached for speed in the comparator
		vector<float>sigmathisv=small_this->calc_radial_dist(ss/2,0,1,4);
		vector<float>sigmatov=small_to->calc_radial_dist(ss/2,0,1,4);
		for (int i=0; i<ss/2; i++) {
			sigmathisv[i]*=sigmathis;
			sigmatov[i]*=sigmato;
		}
		
		// debug out
// 		EMData *x=small_this->do_ift();
// 		x->process_inplace("xform.phaseorigin.tocenter");
// 		x->write_image("dbg.hdf",(sexp-5)*2);
// 		delete x;
// 		x=small_to->do_ift();
// 		x->process_inplace("xform.phaseorigin.tocenter");
// 		x->write_image("dbg.hdf",(sexp-5)*2+1);
// 		delete x;

		// This is a solid estimate for very complete searching, 2.5 is a bit arbitrary
		// make sure the altitude step hits 90 degrees, not absolutely necessary for this, but can't hurt
//		float astep = 89.999/floor(pi/(2.5*2.0*atan(2.0/ss)));
		float astep = 89.999/floor(pi/(1.5*2.0*atan(2.0/ss)));

		// This is drawn from single particle analysis testing, which in that case insures that enough sampling to
		// reasonably fill Fourier space is achieved, but doesn't perfectly apply to SPT
//		float astep = (float)(89.99/ceil(90.0*9.0/(8.0*sqrt((float)(4300.0/ss)))));	// 8 is (3+speed) from SPA with speed=5

		// This insures we make at least one real effort at each level
		for (int i=0; i<nsoln; i++) {
			s_score[i]=1.0e24;	// reset the scores since the different scales will not match
			if (fabs(s_step[i*3+0])<astep/4.0) s_step[i*3+0]*=2.0;
			if (fabs(s_step[i*3+1])<astep/4.0) s_step[i*3+1]*=2.0;
			if (fabs(s_step[i*3+2])<astep/4.0) s_step[i*3+2]*=2.0;
		}
		
		// This is for the first loop, we do a full search in a heavily downsampled space
		if (s_coverage[0]==0.0f) {
			// Genrate points on a sphere in an asymmetric unit
			if (verbose>1) printf("stage 1 - ang step %1.2f\n",astep);
			Dict d;
			d["inc_mirror"] = true;
			d["delta"] = astep;		
			Symmetry3D* sym = Factory<Symmetry3D>::get((string)params.set_default("sym","c1"));
			// We don't generate for phi, since this can produce a very large number of orientations
			vector<Transform> transforms = sym->gen_orientations((string)params.set_default("orientgen","eman"),d);
			if (verbose>0) printf("%d orientations to test (%d)\n",(int)(transforms.size()*(360.0/astep)),transforms.size());
			if (transforms.size()<30) continue; // for very high symmetries we will go up to 32 instead of 24

			// We iterate over all orientations in an asym triangle (alt & az) then deal with phi ourselves
//			for (std::vector<Transform>::iterator t = transforms.begin(); t!=transforms.end(); ++t) {    // iterator form was causing all sorts of problems
			for (unsigned int it=0; it<transforms.size(); it++) {
				Transform t = transforms[it];
				if (verbose>2) {
					printf("  %d/%d \r",it,transforms.size());
					fflush(stdout);
				}
				for (float phi=0; phi<360.0; phi+=astep) {
					Dict aap=t.get_params("eman");
					aap["phi"]=phi;
					aap["tx"]=0;
					aap["ty"]=0;
					aap["tz"]=0;
					t.set_params(aap);
					
					// somewhat strangely, rotations are actually much more expensive than FFTs, so we use a CCF for translation
					EMData *stt=small_this->process("xform",Dict("transform",EMObject(&t),"zerocorners",1));
					EMData *ccf=small_to->calc_ccf(stt);
					IntPoint ml=ccf->calc_max_location_wrap();

					aap["tx"]=(int)ml[0];
					aap["ty"]=(int)ml[1];
					aap["tz"]=(int)ml[2];
					t.set_params(aap);
					delete stt;
					delete ccf;
					stt=small_this->process("xform",Dict("transform",EMObject(&t),"zerocorners",1));	// we have to do 1 slow transform here now that we have the translation

//					float sim=stt->cmp("ccc.tomo.thresh",small_to,Dict("sigmaimg",sigmathis,"sigmawith",sigmato));
					float sim=stt->cmp("ccc.tomo.thresh",small_to);
//					float sim=stt->cmp("fsc.tomo.auto",small_to,Dict("sigmaimg",sigmathisv,"sigmawith",sigmatov));
//					float sim=stt->cmp("fsc.tomo.auto",small_to);

					// We want to make sure our starting points are somewhat separated from each other, so we replace any angles too close to an existing angle
					// If we find an existing 'best' angle within range, then we either replace it or skip
					int worst=-1;
					for (int i=0; i<nsoln; i++) {
						if (s_coverage[i]==0.0) continue;	// hasn't been set yet
						Transform tdif=s_xform[i].inverse();
						tdif=tdif*t;
						float adif=tdif.get_rotation("spin")["omega"];
						if (adif<astep*2.5) {
							worst=i;
//							printf("= %1.3f\n",adif);
						}
					}
					
					// if we weren't close to an existing angle, then we find the lowest current score and use that
					if (worst==-1) {
						// First we find the worst solution in the list of possible best solutions, or the first
						// solution which is currently "empty"
						for (int i=0; i<nsoln; i++) {
							if (s_coverage[i]==0.0) { worst=i; break; }
							if (s_score[i]<s_score[worst]) worst=i;
						}
					}

					// If the current solution is better than the 'worst' of the previous solutions, then we
					// displace it. Note that there is no sorting performed here
					if (sim<s_score[worst]) {
						s_score[worst]=sim;
						s_coverage[worst]=stt->get_attr("fft_overlap");
						s_xform[worst]=t;
						//printf("%f\t%f\t%d\n",s_score[worst],s_coverage[worst],worst);
					}
					delete stt;
				}
			}
			if (verbose>2) printf("\n");


		}
		// Once we have our initial list of best locations, we just refine each possibility individually
		else {
			// We generate a search pattern around each existing solution
			if (verbose>1) printf("stage 2 (%1.2f)\n",astep);
			for (int i=0; i<nsoln; i++) {

				if (verbose>2) {
					printf("  %d\t%d\r",i,nsoln);
					fflush(stdout);
				}
				// We work an axis at a time until we get where we want to be. Somewhat like a simplex
				int changed=1;
				while (changed) {
					changed=0;
					for (int axis=0; axis<3; axis++) {
						if (fabs(s_step[i*3+axis])<astep/4.0) continue;		// skip axes where we already have enough precision on this axis
						Dict upd;
						upd[axname[axis]]=s_step[i*3+axis];
						// when moving az, we move phi in the opposite direction by the same amount since the two are singular at alt=0
						// phi continues to move independently. I believe this should produce a more monotonic energy surface
						if (axis==0) upd[axname[2]]=-s_step[i*3+axis];		

						int r=testort(small_this,small_to,sigmathisv,sigmatov,s_score,s_coverage,s_xform,i,upd);
						
						// If we fail, we reverse direction with a slightly smaller step and try that
						// Whether this fails or not, we move on to the next axis
						if (r) changed=1; 
						else {
							s_step[i*3+axis]*=-0.75;
							upd[axname[axis]]=s_step[i*3+axis];
							r=testort(small_this,small_to,sigmathisv,sigmatov,s_score,s_coverage,s_xform,i,upd);
							if (r) changed=1;
						}
						if (verbose>4) printf("\nX %1.3f\t%1.3f\t%1.3f\t%d\t",s_step[i*3],s_step[i*3+1],s_step[i*3+2],changed);
					}
					if (verbose>3) {
							Dict aap=s_xform[i].get_params("eman");
							printf("\n%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t(%1.3f)",s_step[i*3],s_step[i*3+1],s_step[i*3+2],float(aap["az"]),float(aap["alt"]),float(aap["phi"]),s_score[i]);
					}
					
					if (!changed) {
						for (int j=0; j<3; j++) s_step[i*3+j]*-0.75;
						changed=1;
					}
					if (fabs(s_step[i*3])<astep/4 && fabs(s_step[i*3+1])<astep/4 && fabs(s_step[i*3+2])<astep/4) changed=0;
				}
				
				// Ouch, exhaustive (local) search
// 				for (int daz=-1; daz<=1; daz++) {
// 					for (int dalt=-1; dalt<=1; dalt++) {
// 						for (int dphi=-1; dphi<=1; dphi++) {
// 							Dict upd;
// 							upd["az"]=daz*astep;
// 							upd["alt"]=dalt*astep;
// 							upd["phi"]=dphi*astep;
// 							int r=testort(small_this,small_to,s_score,s_coverage,s_xform,i,upd);
// 						}
// 					}
// 				}
			}
		}
		// lazy earlier in defining s_ vectors, so lazy here too and inefficiently sorting
		// We are sorting inside the outermost loop so we can decrease the number of solutions
		// before we get to the finest precision
		for (unsigned int i=0; i<nsoln-1; i++) {
			for (unsigned int j=i+1; j<nsoln; j++) {
				if (s_score[i]>s_score[j]) {
					float t=s_score[i]; s_score[i]=s_score[j]; s_score[j]=t;
					t=s_coverage[i]; s_coverage[i]=s_coverage[j]; s_coverage[j]=t;
					Transform tt=s_xform[i]; s_xform[i]=s_xform[j]; s_xform[j]=tt;
				}
			}
		}
		
		// At each level of sampling we (potentially) decrease the number of answers we check in detail
		// assuming we are gradually homing in on the best solution
		nsoln/=2;
		if (nsoln<nrsoln) nsoln=nrsoln;

		
		delete small_this;
		delete small_to;
		if (ss==ny) break;
	}
	
	delete base_this;
	delete base_to;


	// initialize results
	vector<Dict> solns;
	for (unsigned int i = 0; i < nrsoln; ++i ) {
		Dict d;
		d["score"] = s_score[i];
		d["coverage"] = s_coverage[i];
		s_xform[i].invert();	// this is because we inverted the order of the input images above to match convention
		d["xform.align3d"] = &s_xform[i];
		solns.push_back(d);
	}

	return solns;
}

// This is just to prevent redundancy. It takes the existing solution vectors as arguments, an a proposed update for
// vector i. It updates the vectors if the proposal makes an improvement, in which case it returns true
bool RT3DTreeAligner::testort(EMData *small_this,EMData *small_to,vector<float> &sigmathisv,vector<float> &sigmatov,vector<float> &s_score, vector<float> &s_coverage,vector<Transform> &s_xform,int i,Dict &upd) const {
	Transform t;
	Dict aap=s_xform[i].get_params("eman");
	aap["tx"]=0;
	aap["ty"]=0;
	aap["tz"]=0;
	for (Dict::const_iterator p=upd.begin(); p!=upd.end(); p++) {
		aap[p->first]=(float)aap[p->first]+(float)p->second;
	}

	t.set_params(aap);

	// rotate in Fourier space then use a CCF to find translation
	EMData *stt=small_this->process("xform",Dict("transform",EMObject(&t),"zerocorners",1));
	EMData *ccf=small_to->calc_ccf(stt);
	IntPoint ml=ccf->calc_max_location_wrap();

	
	aap["tx"]=(int)ml[0];
	aap["ty"]=(int)ml[1];
	aap["tz"]=(int)ml[2];
	t.set_params(aap);
	EMData *st2=small_this->process("xform",Dict("transform",EMObject(&t),"zerocorners",1));	// we have to do 1 slow transform here now that we have the translation
	
	float sim=st2->cmp("fsc.tomo.auto",small_to,Dict("sigmaimg",sigmathisv,"sigmawith",sigmatov));
//	float sim=st2->cmp("ccc.tomo.thresh",small_to,Dict("sigmaimg",sigmathis,"sigmawith",sigmato));
// 	printf("\nTESTORT %6.1f  %6.1f  %6.1f\t%4d %4d %4d\t%1.5g\t%1.5g %d (%d)",
// 		float(aap["az"]),float(aap["alt"]),float(aap["phi"]),int(aap["tx"]),int(aap["ty"]),int(aap["tz"]),sim,s_score[i],int(sim<s_score[i]),ccf->get_ysize());

	delete ccf;
	// If the score is better than before, we update this particular best value
	if (sim<s_score[i]) {
		s_score[i]=sim;
		s_coverage[i]=st2->get_attr("fft_overlap");
		s_xform[i]=t;
//  		if (sim<-.67 && stt->get_ysize()==160) {
//  			stt->write_image("dbug2.hdf",0);
//  			st2->write_image("dbug2.hdf",1);
//  			small_to->write_image("dbug2.hdf",2);
//  			printf("$$$$\n");
// 			t.printme();
//  		}
		delete stt;
		delete st2;
		return true;
	}
	delete stt;
	delete st2;
	return false;
}


EMData* RT3DSphereAligner::align(EMData * this_img, EMData *to, const string & cmp_name, const Dict& cmp_params) const
{

	vector<Dict> alis = xform_align_nbest(this_img,to,1,cmp_name,cmp_params);

	Dict t;
	Transform* tr = (Transform*) alis[0]["xform.align3d"];
	t["transform"] = tr;
	EMData* soln = this_img->process("xform",t);
	soln->set_attr("xform.align3d",tr);
	delete tr; tr = 0;

	return soln;

}

vector<Dict> RT3DSphereAligner::xform_align_nbest(EMData * this_img, EMData * to, const unsigned int nsoln, const string & cmp_name, const Dict& cmp_params) const {

	if ( this_img->get_ndim() != 3 || to->get_ndim() != 3 ) {
		throw ImageDimensionException("This aligner only works for 3D images");
	}

	int searchx = 0;
	int searchy = 0;
	int searchz = 0;

	bool dotrans = params.set_default("dotrans",1);
	if (params.has_key("search")) {
		vector<string> check;
		check.push_back("searchx");
		check.push_back("searchy");
		check.push_back("searchz");
		for(vector<string>::const_iterator cit = check.begin(); cit != check.end(); ++cit) {
			if (params.has_key(*cit)) throw InvalidParameterException("The search parameter is mutually exclusive of the searchx, searchy, and searchz parameters");
		}
		int search  = params["search"];
		searchx = search;
		searchy = search;
		searchz = search;
	} else {
		searchx = params.set_default("searchx",3);
		searchy = params.set_default("searchy",3);
		searchz = params.set_default("searchz",3);
	}

	Transform* initxform;
	if (params.has_key("initxform") ) {
		// Unlike the 2d refine aligner, this class doesn't require the starting transform's
		// parameters to form the starting guess. Instead the Transform itself
		// is perturbed carefully (using quaternion rotation) to overcome problems that arise
		// when you use orthogonally-based Euler angles
		initxform = params["initxform"];
	}else {
		initxform = new Transform(); // is the identity
	}

	float lphi = params.set_default("phi0",0.0f);
	float uphi = params.set_default("phi1",360.0f);
	float dphi = params.set_default("dphi",10.f);
	float threshold = params.set_default("threshold",0.f);
	if (threshold < 0.0f) throw InvalidParameterException("The threshold parameter must be greater than or equal to zero");
	bool verbose = params.set_default("verbose",false);

	//in case we are aligning tomos
	Dict altered_cmp_params(cmp_params);
	if (cmp_name == "ccc.tomo") {
		altered_cmp_params.set_default("searchx", searchx);
		altered_cmp_params.set_default("searchy", searchy);
		altered_cmp_params.set_default("searchz", searchz);
		altered_cmp_params.set_default("norm", true);
	}

	vector<Dict> solns;
	if (nsoln == 0) return solns; // What was the user thinking?
	for (unsigned int i = 0; i < nsoln; ++i ) {
		Dict d;
		d["score"] = 1.e24;
		Transform t; // identity by default
		d["xform.align3d"] = &t; // deep copy is going on here
		solns.push_back(d);
	}

	Dict d;
	d["inc_mirror"] = true; // This should probably always be true for 3D case. If it ever changes we might have to make inc_mirror a parameter
	if ( params.has_key("delta") && params.has_key("n") ) {
		throw InvalidParameterException("The delta and n parameters are mutually exclusive in the RT3DSphereAligner aligner");
	} else if ( params.has_key("n") ) {
		d["n"] = params["n"];
	} else {
		// If they didn't specify either then grab the default delta - if they did supply delta we're still safe doing this
		d["delta"] = params.set_default("delta",10.f);
	}

	if ((string)params.set_default("orientgen","eman")=="eman") d["perturb"]=0;
	Symmetry3D* sym = Factory<Symmetry3D>::get((string)params.set_default("sym","c1"));
	vector<Transform> transforms = sym->gen_orientations((string)params.set_default("orientgen","eman"),d);

	bool tomography = (cmp_name == "ccc.tomo") ? 1 : 0;

	//precompute fixed FT, saves a LOT of time!!!
	EMData * this_imgfft = 0;
	if(dotrans || tomography){
		this_imgfft = this_img->do_fft();
	}

#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		cout << "Using CUDA for 3D alignment" << endl;
		if(!to->getcudarodata()) to->copy_to_cudaro(); // Safer call
		if(!this_img->getcudarwdata()) this_img->copy_to_cuda();
		if(this_imgfft) this_imgfft->copy_to_cuda();
	}
#endif

	Transform trans = Transform();
	Cmp* c = Factory <Cmp>::get(cmp_name, cmp_params);

	bool use_cpu = true;
	for(vector<Transform>::const_iterator trans_it = transforms.begin(); trans_it != transforms.end(); trans_it++) {
		Dict params = trans_it->get_params("eman");

		if (verbose) {
			float alt = params["alt"];
			float az = params["az"];
			cout << "Trying angle alt: " << alt << " az: " << az << endl;
		}

		for( float phi = lphi; phi < uphi; phi += dphi ) {
			params["phi"] = phi;
			Transform t(params);
			t = t*(*initxform);

			EMData* transformed;
			transformed = to->process("xform",Dict("transform",&t));

			//need to do things a bit diffrent if we want to compare two tomos
			float best_score = 0.0f;
			// Dotrans is effectievly ignored for tomography
			if(dotrans || tomography){
				EMData* ccf = transformed->calc_ccf(this_imgfft);
#ifdef EMAN2_USING_CUDA
				if(EMData::usecuda == 1){
					// I use the following code rather than ccc.tomo to avoid doing two CCCs
					use_cpu = false;
					CudaPeakInfo* data = calc_max_location_wrap_cuda(ccf->getcudarwdata(), ccf->get_xsize(), ccf->get_ysize(), ccf->get_zsize(), searchx, searchy, searchz);
					trans.set_trans((float)-data->px, (float)-data->py, (float)-data->pz);
					t = trans*t;	//composite transform to reflect the fact that we have done a rotation first and THEN a transformation
					if (tomography) {
						float2 stats = get_stats_cuda(ccf->getcudarwdata(), ccf->get_xsize(), ccf->get_ysize(), ccf->get_zsize());
						best_score = -(data->peak - stats.x)/sqrt(stats.y); // Normalize, this is better than calling the norm processor since we only need to normalize one point
					} else {
						best_score = -data->peak;
					}
					delete data;
				}
#endif
				if(use_cpu){
					// I use the following code rather than ccc.tomo to avoid doing two CCCs
					if(tomography) ccf->process_inplace("normalize");
					IntPoint point = ccf->calc_max_location_wrap(searchx,searchy,searchz);
					trans.set_trans((float)-point[0], (float)-point[1], (float)-point[2]);
					t = trans*t;	//composite transform to reflect the fact that we have done a rotation first and THEN a transformation
					best_score = -ccf->get_value_at_wrap(point[0], point[1], point[2]);
				}
				delete ccf; ccf =0;
				delete transformed; transformed = 0;// this is to stop a mem leak
			}

			if(!tomography){
				if(!transformed) transformed = to->process("xform",Dict("transform",&t));
				best_score = c->cmp(this_img,transformed);
				delete transformed; transformed = 0;
			}

			unsigned int j = 0;
			//cout << "alt " <<float(t.get_rotation("eman").get("alt")) << " az " << float(t.get_rotation("eman").get("az")) << " phi " << float(t.get_rotation("eman").get("phi")) << endl;
			for ( vector<Dict>::iterator it = solns.begin(); it != solns.end(); ++it, ++j ) {
				if ( (float)(*it)["score"] > best_score ) { // Note greater than - EMAN2 preferes minimums as a matter of policy
					vector<Dict>::reverse_iterator rit = solns.rbegin();
					copy(rit+1,solns.rend()-j,rit);
					Dict& d = (*it);
					d["score"] = best_score;
					t.invert(); //We actually moved the ref onto the moving, so we need to invert to do the opposite(this is done b/c the ref is aligned to the sym axis, whereas the mvoing is not)
					d["xform.align3d"] = &t; // deep copy is going on here
					break;
				}
			}

		}
	}

	if(this_imgfft) {delete this_imgfft; this_imgfft = 0;}
	if(sym!=0) delete sym;
	if (c != 0) delete c;

	return solns;

}

//Could refactor the code here......(But not really woth it)
EMData* RT3DSymmetryAligner::align(EMData * this_img, EMData *to, const string & cmp_name, const Dict& cmp_params) const
{

	vector<Dict> alis = xform_align_nbest(this_img,to,1,cmp_name,cmp_params);

	Transform* tr = (Transform*) alis[0]["xform.align3d"];
	EMData* soln = this_img->process("xform",Dict("transform",tr));
	soln->set_attr("xform.align3d",tr);
	delete tr; tr = 0;

	return soln;

}

vector<Dict> RT3DSymmetryAligner::xform_align_nbest(EMData * this_img, EMData * to, const unsigned int nsoln, const string & cmp_name, const Dict& cmp_params) const
{

	bool verbose = params.set_default("verbose",false);
	Transform* ixform;
	if (params.has_key("transform") ) {
		ixform = params["transform"];
	}else{
		ixform = new Transform(); // is the identity
	}

	//Initialize a soln dict
	vector<Dict> solns;
	if (nsoln == 0) return solns; // What was the user thinking?
	for (unsigned int i = 0; i < nsoln; ++i ) {
		Dict d;
		d["score"] = 1.e24;
		Transform t; // identity by default
		d["xform.align3d"] = &t; // deep copy is going on here
		solns.push_back(d);
	}

	#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1) {
		cout << "Using CUDA for 3D sym alignment" << endl;
		if(!this_img->getcudarwdata()) this_img->copy_to_cudaro();
		if(!to->getcudarwdata()) to->copy_to_cuda();
	}
	#endif

	//Generate symmetry related orientations
	vector<Transform> syms = Symmetry3D::get_symmetries((string)params.set_default("sym","icos"));
	Cmp* c = Factory <Cmp>::get(cmp_name, cmp_params);

	float score = 0.0f;
	for ( vector<Transform>::const_iterator symit = syms.begin(); symit != syms.end(); ++symit ) {
		//Here move to sym position and compute the score
		Transform sympos = (*symit)*(*ixform);
		EMData* transformed = this_img->process("xform",Dict("transform", &sympos));
		score = c->cmp(transformed,to);
		delete transformed; transformed = 0;

		if (verbose) {
			Dict rots = sympos.get_rotation("eman");
			cout <<"Score is: " << score << " az " << float(rots["az"]) << " alt " << float(rots["alt"]) << " phi " << float(rots["phi"]) << endl;
		}

		unsigned int j = 0;
		for ( vector<Dict>::iterator it = solns.begin(); it != solns.end(); ++it, ++j ) {
			if ( (float)(*it)["score"] > score ) { // Note greater than - EMAN2 preferes minimums as a matter of policy
				vector<Dict>::reverse_iterator rit = solns.rbegin();
				copy(rit+1,solns.rend()-j,rit);
				Dict& d = (*it);
				d["score"] = score;
				d["xform.align3d"] = &sympos; // deep copy is going on here
				break;
			}
		}
	}

	if (c != 0) delete c;

	return solns;
}

namespace {
float frm_2d_Align(EMData *this_img, EMData *to, float *frm2dhhat, EMData* selfpcsfft, int p_max_input,int rsize, float &com_this_x, float &com_this_y, float &com_with_x, float &com_with_y,const string & cmp_name, const Dict& cmp_params)
{
	int size=rsize;
	float dx,dy;
	int bw=size/2;
	int MAXR=this_img->get_ysize()/2;
	//int MAXR=size;
	unsigned long tsize=2*size;
	unsigned long ind1=0, ind2=0, ind3=0, ind4=0, ind41=0;
	unsigned long index0=0;
	int i=0, j=0, m=0, n=0, r=0;
	int loop_rho=0, rho_best=0;

 	float* gnr2   = new float[size*2];
 	float* maxcor = new float[size+1];                  // MAXR need change

 	int p_max=p_max_input;
	float* result = new float[5*(p_max+1)];
	float* cr=new float[size*(bw+1)];
	float* ci=new float[size*(bw+1)];
	EMData *data_in=new EMData;
	data_in->set_complex(true);
	data_in->set_fftodd(false);
	data_in->set_ri(true);
	data_in->set_size(size+2,size,1);
	float *in=data_in->get_data();

	float *self_sampl_fft = selfpcsfft->get_data(); // ming f(r)

	float maxcor_sofar=0.0f;
	int p=0;

	for(p=0; p<=p_max; ++p){
		ind1=p*size*bw;
		for (i=0;i<size;++i)
			for (j=0;j<bw+1;++j){
				cr[i*(bw+1)+j]=0.0;
				ci[i*(bw+1)+j]=0.0;
			}
    	for(n=0;n<bw;++n){                                // loop for n
    		ind2=(ind1+n);
      		index0=n*(bw+1);
			for(r=0;r<=MAXR;++r) {
      			ind3=(ind2+r*bw)*size;
      			for(m=0;m<size;m++){              // take back hat{h(n,r,p)}(m)
      				ind4=(ind3+m)*2;
				    ind41=ind4+1;
				    gnr2[2*m]=frm2dhhat[ind4];
				    gnr2[2*m+1]=frm2dhhat[ind41];
				}
      			for(m=0;m<bw;++m){
					float tempr=self_sampl_fft[2*m+r*(size+2)]*r;
      				float tempi=self_sampl_fft[2*m+1+r*(size+2)]*r;
      				float gnr2_r=gnr2[2*m];
      				float gnr2_i=gnr2[2*m+1];
      				cr[n*(bw+1)+m]+=gnr2_r*tempr+gnr2_i*tempi;
					ci[n*(bw+1)+m]+=gnr2_i*tempr-gnr2_r*tempi;
					if(n!=0){					// m,-n
      				   	if(m!= 0){
				      		int ssize=tsize-2*m;	// ssize = 2*size-2m
				      		int ssize1=ssize+1;
				      		float gnr2_r=gnr2[ssize];
				      		float gnr2_i=gnr2[ssize1];
					   		cr[(size-n)*(bw+1)+m]+=gnr2_r*tempr-gnr2_i*tempi;
				      		ci[(size-n)*(bw+1)+m]-=gnr2_i*tempr+gnr2_r*tempi;
				      	}
				   		else{
				   			cr[(size-n)*(bw+1)+m]+=*(gnr2)*tempr-*(gnr2+1)*tempi;
				   			ci[(size-n)*(bw+1)+m]-=*(gnr2+1)*tempr+*(gnr2)*tempi;
				   		}
			       	}
				}
			}
        }
    	for (int cii=0; cii<size*(bw+1);++cii){
    			in[2*cii]=cr[cii];
    			in[2*cii+1]=ci[cii];
    			//printf("cii=%d,in[2i+1]=%f\n",cii, cr[cii]);
    	}

    	EMData *data_out;
		data_out=data_in->do_ift();
		float *c=data_out->get_data();
		float tempr=0.0f, corre_fcs=999.0f;

  	    int n_best=0, m_best=0;
        float temp=-100.0f;
  		for(n=0;n<size;++n){// move Tri_2D to Tri = c(phi,phi';rho)
  		  	for(m=0;m<size;++m){
				temp=c[n*size+m];
				if(temp>tempr) {
					tempr=temp;
					n_best=n;
					m_best=m;
				}
  		   	}
  		}
  		delete data_out;

  		float corre,Phi2,Phi,Tx,Ty,Vx, Vy;

  		//for (n_best=0;n_best<bw;n_best++)
  		  //  for (m_best=0;m_best<2*bw;m_best++){
  		//n_best=0;
  		//m_best=70;
  		Phi2=n_best*M_PI/bw;  // ming this is reference image rotation angle
  		Phi=m_best*M_PI/bw;   // ming this is particle image rotation angle
  		Vx=p*cos(Phi);//should use the angle of the centered one
  		Vy=-p*sin(Phi);
  		Tx=Vx+(floor(com_this_x+0.5f)-floor(com_with_x+0.5f));
  		Ty=Vy+(floor(com_this_y+0.5f)-floor(com_with_y+0.5f));

  		dx=-Tx;	// the Rota & Trans value (Tx,Ty, ang_keep) are for the projection image,
  		dy=-Ty;	// need to convert to raw image

  		EMData *this_tmp=this_img->copy();//ming change to to
		this_tmp->rotate(-(Phi2-Phi)*180.0f,0.0f,0.0f);
		this_tmp->translate(dx,dy,0.0);

		corre=this_tmp->cmp(cmp_name,to,cmp_params);
		//printf("corre=%f\n",corre);
		delete this_tmp;
		if(corre<=corre_fcs) { //ming, cmp use smaller value stands for more similarity
			corre_fcs=corre;
			result[0+5*p] = float(p);	// rho
			result[1+5*p] = corre;		// correlation_fcs
			result[2+5*p] = (Phi2-Phi)*180.0f;	// rotation angle
			result[3+5*p] = Tx;		// Translation_x
			result[4+5*p] = Ty;		// Translation_y
		}
		maxcor[p]=corre_fcs;               		//  maximum correlation for current rho
		if(corre_fcs<maxcor_sofar) {
			maxcor_sofar=corre_fcs;   		// max correlation up to current rho
		    rho_best=p;				// the rho value with maxinum correlation value
		}
		if(p>=4){
			if(maxcor[p] < maxcor[p-1] && maxcor[p-1] < maxcor[p-2]&& maxcor[p-2] < maxcor[p-3] && maxcor[p-3] < maxcor[p-4]){
				loop_rho=1;
				break; //exit p loop
			}
		}
	} // end for p
	//}//test my method
	if(loop_rho == 1)
		p=p+1;
	int rb5=5*rho_best;
	float fsc      = result[1+rb5];
	float ang_keep = result[2+rb5];
	float Tx       = result[3+rb5];
	float Ty       = result[4+rb5];
	delete[] gnr2;
	delete[] maxcor;
	delete[] result;
	delete cr;
	cr=0;
	delete ci;
	ci=0;
	delete data_in; // ming add
	dx = -Tx;		// the Rota & Trans value (Tx,Ty, ang_keep) are for the projection image,
	dy = -Ty;		// need to convert to raw image
	this_img->rotate(-ang_keep,0,0); // ming change this to this_img??
	this_img->translate(dx,dy,0.0); // ming change this to this_img


	Transform  tsoln(Dict("type","2d","alpha",ang_keep));
	tsoln.set_trans(dx,dy);
	this_img->set_attr("xform.align2d",&tsoln);
#ifdef DEBUG
	float fsc_best=this_img->cmp(cmp_name,to,cmp_params);
	printf("rho_best=%d fsc=%f fsc_best=%f dx=%f dy=%f ang_keep=%f com_withx=%f com_selfx=%f com_selfy=%f\n",rho_best,fsc,fsc_best,dx,dy,ang_keep,com_with_x,com_this_x,com_this_y);
#endif
	return fsc;     // return the fsc coefficients
} // FRM2D aligner sub_class
} // end namespace


EMData *FRM2DAligner::align(EMData * this_img, EMData * to,
			const string & cmp_name, const Dict& cmp_params) const
{
	if (!this_img) {
		return 0;
	}
	if (to && !EMUtil::is_same_size(this_img, to))
		throw ImageDimensionException("Images must be the same size to perform translational alignment");

	int nx=this_img->get_xsize();
	int ny=this_img->get_ysize();
	int size =(int)floor(M_PI*ny/4.0);
	size =Util::calc_best_fft_size(size);//ming   bestfftsize(size);
	int MAXR=ny/2;
	//int MAXR=size;
	EMData *this_temp=this_img->copy(); // ming change avg to to
	FloatPoint com_test,com_test1;
	com_test=this_temp->calc_center_of_mass();//ming add
	float com_this_x=com_test[0];
	float com_this_y=com_test[1];
	delete this_temp;


	EMData *that_temp=to->copy();
	com_test1=that_temp->calc_center_of_mass();
	float com_with_x=com_test1[0];
	float com_with_y=com_test1[1];
	delete that_temp;

	EMData *avg_frm=to->copy();
	float dx,dy;
	//float dx=-(com_with_x-nx/2); //ming
	//float dy=-(com_with_y-ny/2); //ming
	//avg_frm->translate(dx,dy,0.0);
	EMData *withpcs=avg_frm->unwrap_largerR(0,MAXR,size,float(MAXR)); // ming, something wrong inside this subroutine
	//EMData *withpcs=avg_frm->unwrap(-1,-1,-1,0,0,1);
	EMData *withpcsfft=withpcs->oneDfftPolar(size, float(MAXR), float(MAXR));

	float *sampl_fft=withpcsfft->get_data(); //
	delete avg_frm;
	delete withpcs;

	int bw=size/2;
	unsigned long ind1=0, ind2=0, ind3=0, ind4=0, ind41=0;
	float pi2=2.0*M_PI, r2;

	EMData *data_in=new EMData;
	data_in->set_complex(true);
	data_in->set_ri(1);
	data_in->set_size(2*size,1,1);
	float * comp_in=data_in->get_data();

	int p_max=3;
	float *frm2dhhat=0;

	if( (frm2dhhat=(float *)malloc((p_max+1)*(size+2)*bw*size*2* sizeof(float)))==NULL){
		cout <<"Error in allocating memory 13. \n";
		exit(1);
	}
	//printf("p_max=%d\n",p_max);
	float *sb=0, *cb=0;		// sin(beta) and cos(beta) for get h_hat, 300>size
	if((sb=new float[size])==NULL||(cb=new float[size])==NULL) {
		cout <<"can't allocate more memory, terminating. \n";
		exit(1);
	}
	for(int i=0;i<size;++i) {        // beta sampling, to calculate beta' and r'
		float beta=i*M_PI/bw;
	   	sb[i]=sin(beta);
	   	cb[i]=cos(beta);
	}

	for(int p=0; p<=p_max; ++p){
		ind1=p*size*bw;
    	float pp2=(float)(p*p);
   		for(int n=0;n<bw;++n){         /* loop for n */
    		ind2=ind1+n;
      		for(int r=0;r<=MAXR;++r) {
				ind3=(ind2+r*bw)*size;
        		float rr2=(float)(r*r);
				float rp2=(float)(r*p);
       			for(int i=0;i<size;++i){                            // beta sampling, to get beta' and r'
       				r2=std::sqrt((float)(rr2+pp2-2.0*rp2*cb[i]));   // r2->r'
       		 		int r1=(int)floor(r2+0.5f);                        // for computing gn(r')
       				if(r1>MAXR){
       					comp_in[2*i]=0.0f;
       					comp_in[2*i+1]=0.0f;
       				}
       				else{
       					float gn_r=sampl_fft[2*n+r1*(size+2)];           // real part of gn(r')
       					float gn_i=sampl_fft[2*n+1+r1*(size+2)];           // imaginary part of gn(r')
						float sb2, cb2, cn, sn;
						if(n!=0){
							if(r2 != 0.0){
								sb2=r*sb[i]/r2;
								cb2=(r*cb[i]-p)/r2;
							}
        					else{
								sb2=0.0;
								cb2=1.0;
							}
        					if(sb2>1.0) sb2= 1.0f;
        					if(sb2<-1.0)sb2=-1.0f;
        					if(cb2>1.0) cb2= 1.0f;
        					if(cb2<-1.0)cb2=-1.0f;
        					float beta2=atan2(sb2,cb2);
        					if(beta2<0.0) beta2+=pi2;
        					float nb2=n*beta2;
        					cn=cos(nb2);
							sn=sin(nb2);
						}
        				else{
							cn=1.0f; sn=0.0f;
						}
						comp_in[2*i]=cn*gn_r-sn*gn_i;
						comp_in[2*i+1]=-(cn*gn_i+sn*gn_r);
        			}
        		}
       			EMData *data_out;
        		data_out=data_in->do_fft();
        		float * comp_out=data_out->get_data();
        		for(int m=0;m<size;m++){                                     // store hat{h(n,r,p)}(m)
					ind4=(ind3+m)*2;
					ind41=ind4+1;
					frm2dhhat[ind4]=comp_out[2*m];
					frm2dhhat[ind41]=comp_out[2*m+1];
				}
        		delete data_out;
			}
		}
	}

	delete[] sb;
	delete[] cb;
	delete data_in;
	delete withpcsfft;

	float dot_frm0=0.0f, dot_frm1=0.0f;
	EMData *da_nFlip=0, *da_yFlip=0, *dr_frm=0;
	//dr_frm=this_img->copy();
	for (int iFlip=0;iFlip<=1;++iFlip){
		if (iFlip==0){dr_frm=this_img->copy(); 	da_nFlip=this_img->copy();}
		else {dr_frm=this_img->copy(); da_yFlip=this_img->copy();}
		if(iFlip==1) {com_this_x=nx-com_this_x; } //ming   // image mirror about Y axis, so y keeps the same

		dx=-(com_this_x-nx/2); //ming
		dy=-(com_this_y-ny/2); //ming
		dr_frm->translate(dx,dy,0.0); // this
		EMData *selfpcs = dr_frm->unwrap_largerR(0,MAXR,size, (float)MAXR);
		//EMData *selfpcs=dr_frm->unwrap(-1,-1,-1,0,0,1);
		EMData *selfpcsfft = selfpcs->oneDfftPolar(size, (float)MAXR, (float)MAXR);
		delete selfpcs;
		delete dr_frm;
		if(iFlip==0)
			dot_frm0=frm_2d_Align(da_nFlip,to, frm2dhhat, selfpcsfft, p_max, size, com_this_x, com_this_y, com_with_x, com_with_y,cmp_name,cmp_params);
		else
			dot_frm1=frm_2d_Align(da_yFlip,to, frm2dhhat, selfpcsfft, p_max, size, com_this_x, com_this_y, com_with_x, com_with_y,cmp_name,cmp_params);
		delete selfpcsfft;
	}

	delete[] frm2dhhat;
	if(dot_frm0 <=dot_frm1) {
#ifdef DEBUG
		printf("best_corre=%f, no flip\n",dot_frm0);
#endif
		delete da_yFlip;
		return da_nFlip;
	}
	else {
#ifdef DEBUG
		printf("best_corre=%f, flipped\n",dot_frm1);
#endif
		delete da_nFlip;
		return da_yFlip;
	}
}

#ifdef SPARX_USING_CUDA

CUDA_Aligner::CUDA_Aligner(int id) {
	image_stack = NULL;
	image_stack_filtered = NULL;
	ccf = NULL;
	if (id != -1) cudasetup(id);
}

void CUDA_Aligner::finish() {
	if (image_stack) free(image_stack);
	if (image_stack_filtered) free(image_stack_filtered);
	if (ccf) free(ccf);
	image_stack = NULL;
	image_stack_filtered = NULL;
	ccf = NULL;
}

void CUDA_Aligner::setup(int nima, int nx, int ny, int ring_length, int nring, int ou, float step, int kx, int ky, bool ctf) {

	NIMA = nima;
	NX = nx;
	NY = ny;
	RING_LENGTH = ring_length;
        NRING = nring;
	STEP = step;
	KX = kx;
	KY = ky;
	OU = ou;
	CTF = ctf;

	image_stack = (float *)malloc(NIMA*NX*NY*sizeof(float));
	if (CTF == 1) image_stack_filtered = (float *)malloc(NIMA*NX*NY*sizeof(float));
	ccf = (float *)malloc(2*(2*KX+1)*(2*KY+1)*NIMA*(RING_LENGTH+2)*sizeof(float));
}

void CUDA_Aligner::insert_image(EMData *image, int num) {

	int base_address = num*NX*NY;

	for (int y=0; y<NY; y++)
		for (int x=0; x<NX; x++)
			image_stack[base_address+y*NX+x] = (*image)(x, y);
}

void CUDA_Aligner::filter_stack(vector<float> ctf_params) {

	float *params;

	params = (float *)malloc(NIMA*6*sizeof(float));

	for (int i=0; i<NIMA*6; i++) params[i] = ctf_params[i];

	filter_image(image_stack, image_stack_filtered, NIMA, NX, NY, params);

	free(params);
}

void CUDA_Aligner::sum_oe(vector<float> ctf_params, vector<float> ali_params, EMData *ave1, EMData *ave2) {

	float *ctf_p, *ali_p, *av1, *av2;

	ctf_p = (float *)malloc(NIMA*6*sizeof(float));
	ali_p = (float *)malloc(NIMA*4*sizeof(float));

	if (CTF == 1) {
		for (int i=0; i<NIMA*6; i++)  ctf_p[i] = ctf_params[i];
	}
	for (int i=0; i<NIMA*4; i++)   ali_p[i] = ali_params[i];

	av1 = ave1->get_data();
	av2 = ave2->get_data();

	rot_filt_sum(image_stack, NIMA, NX, NY, CTF, ctf_p, ali_p, av1, av2);

	free(ctf_p);
	free(ali_p);
}

vector<float> CUDA_Aligner::alignment_2d(EMData *ref_image_em, vector<float> sx_list, vector<float> sy_list, int silent) {

	float *ref_image, max_ccf;
	int base_address, ccf_offset;
	float ts, tm;
	float ang, sx = 0, sy = 0, mirror, co, so, sxs, sys;
	float *sx2, *sy2;
	vector<float> align_result;

	sx2 = (float *)malloc(NIMA*sizeof(float));
	sy2 = (float *)malloc(NIMA*sizeof(float));

	ref_image = ref_image_em->get_data();

	for (int i=0; i<NIMA; i++) {
		sx2[i] = sx_list[i];
		sy2[i] = sy_list[i];
	}

	if (CTF == 1) {
		calculate_ccf(image_stack_filtered, ref_image, ccf, NIMA, NX, NY, RING_LENGTH, NRING, OU, STEP, KX, KY, sx2, sy2, silent);
	} else {
		calculate_ccf(image_stack, ref_image, ccf, NIMA, NX, NY, RING_LENGTH, NRING, OU, STEP, KX, KY, sx2, sy2, silent);
	}

	ccf_offset = NIMA*(RING_LENGTH+2)*(2*KX+1)*(2*KY+1);

	for (int im=0; im<NIMA; im++) {
		max_ccf = -1.0e22;
		for (int kx=-KX; kx<=KX; kx++) {
			for (int ky=-KY; ky<=KY; ky++) {
				base_address = (((ky+KY)*(2*KX+1)+(kx+KX))*NIMA+im)*(RING_LENGTH+2);
				for (int l=0; l<RING_LENGTH; l++) {
					ts = ccf[base_address+l];
					tm = ccf[base_address+l+ccf_offset];
					if (ts > max_ccf) {
						ang = float(l)/RING_LENGTH*360.0;
						sx = -kx*STEP;
						sy = -ky*STEP;
						mirror = 0;
						max_ccf = ts;
					}
					if (tm > max_ccf) {
						ang = float(l)/RING_LENGTH*360.0;
						sx = -kx*STEP;
						sy = -ky*STEP;
						mirror = 1;
						max_ccf = tm;
					}
				}
			}
		}
		co =  cos(ang*M_PI/180.0);
		so = -sin(ang*M_PI/180.0);
		sxs = sx*co - sy*so;
		sys = sx*so + sy*co;

		align_result.push_back(ang);
		align_result.push_back(sxs);
		align_result.push_back(sys);
		align_result.push_back(mirror);
	}

	free(sx2);
	free(sy2);

	return align_result;
}


vector<float> CUDA_Aligner::ali2d_single_iter(EMData *ref_image_em, vector<float> ali_params, float csx, float csy, int silent, float delta) {

	float *ref_image, max_ccf;
	int base_address, ccf_offset;
	float ts, tm;
	float ang = 0.0, sx = 0.0, sy = 0.0, co, so, sxs, sys;
	int mirror;
	float *sx2, *sy2;
	vector<float> align_result;

	sx2 = (float *)malloc(NIMA*sizeof(float));
	sy2 = (float *)malloc(NIMA*sizeof(float));

	ref_image = ref_image_em->get_data();

	for (int i=0; i<NIMA; i++) {
		ang = ali_params[i*4]/180.0*M_PI;
		sx = (ali_params[i*4+3] < 0.5)?(ali_params[i*4+1]-csx):(ali_params[i*4+1]+csx);
		sy = ali_params[i*4+2]-csy;
		co = cos(ang);
		so = sin(ang);
		sx2[i] = -(sx*co-sy*so);
		sy2[i] = -(sx*so+sy*co);
	}

	if (CTF == 1) {
		calculate_ccf(image_stack_filtered, ref_image, ccf, NIMA, NX, NY, RING_LENGTH, NRING, OU, STEP, KX, KY, sx2, sy2, silent);
	} else {
		calculate_ccf(image_stack, ref_image, ccf, NIMA, NX, NY, RING_LENGTH, NRING, OU, STEP, KX, KY, sx2, sy2, silent);
	}

	ccf_offset = NIMA*(RING_LENGTH+2)*(2*KX+1)*(2*KY+1);

	float sx_sum = 0.0f;
	float sy_sum = 0.0f;

	int dl;
	dl = static_cast<int>(delta/360.0*RING_LENGTH);
	if (dl<1) { dl = 1; }

	for (int im=0; im<NIMA; im++) {
		max_ccf = -1.0e22;
		for (int kx=-KX; kx<=KX; kx++) {
			for (int ky=-KY; ky<=KY; ky++) {
				base_address = (((ky+KY)*(2*KX+1)+(kx+KX))*NIMA+im)*(RING_LENGTH+2);
				for (int l=0; l<RING_LENGTH; l+=dl) {
					ts = ccf[base_address+l];
					tm = ccf[base_address+l+ccf_offset];
					if (ts > max_ccf) {
						ang = float(l)/RING_LENGTH*360.0;
						sx = -kx*STEP;
						sy = -ky*STEP;
						mirror = 0;
						max_ccf = ts;
					}
					if (tm > max_ccf) {
						ang = float(l)/RING_LENGTH*360.0;
						sx = -kx*STEP;
						sy = -ky*STEP;
						mirror = 1;
						max_ccf = tm;
					}
				}
			}
		}
		co =  cos(ang*M_PI/180.0);
		so = -sin(ang*M_PI/180.0);

		sxs = (sx-sx2[im])*co-(sy-sy2[im])*so;
		sys = (sx-sx2[im])*so+(sy-sy2[im])*co;

		//if (sxs*sxs+sys*sys >= 7*7) { sxs=0; sys=0; }

		align_result.push_back(ang);
		align_result.push_back(sxs);
		align_result.push_back(sys);
		align_result.push_back(mirror);

		if (mirror == 0)  { sx_sum += sxs; }  else { sx_sum -= sxs; }
		sy_sum += sys;
	}

	align_result.push_back(sx_sum);
	align_result.push_back(sy_sum);

	free(sx2);
	free(sy2);

	return align_result;
}


CUDA_multiref_aligner::CUDA_multiref_aligner(int id) {
	image_stack = NULL;
	ref_image_stack = NULL;
	ref_image_stack_filtered = NULL;
	ccf = NULL;
	ctf_params = NULL;
	ali_params = NULL;
	cudasetup(id);
}


void CUDA_multiref_aligner::finish() {
	if (image_stack) free(image_stack);
	if (ref_image_stack) free(ref_image_stack);
	if (ref_image_stack_filtered) free(ref_image_stack_filtered);
	if (ccf) free(ccf);
	if (ctf_params) free(ctf_params);
	if (ali_params) free(ali_params);
	image_stack = NULL;
	ref_image_stack = NULL;
	ref_image_stack_filtered = NULL;
	ccf = NULL;
	ctf_params = NULL;
	ali_params = NULL;
}

void CUDA_multiref_aligner::setup(int nima, int nref, int nx, int ny, int ring_length, int nring, int ou, float step, int kx, int ky, bool ctf) {

	NIMA = nima;
	NREF = nref;
	NX = nx;
	NY = ny;
	RING_LENGTH = ring_length;
        NRING = nring;
	STEP = step;
	KX = kx;
	KY = ky;
	OU = ou;
	CTF = ctf;
	// This number can be increased according to the GPU memory. But my tests has shown the speedup
	// is limited (~5%) even if I increased the size 10 times, so it's better to be on the safe side.
	MAX_IMAGE_BATCH = 10;

	image_stack = (float *)malloc(NIMA*NX*NY*sizeof(float));
	ref_image_stack = (float *)malloc(NREF*NX*NY*sizeof(float));
	if (CTF == 1) ref_image_stack_filtered = (float *)malloc(NREF*NX*NY*sizeof(float));
	ccf = (float *)malloc(2*(2*KX+1)*(2*KY+1)*NREF*(RING_LENGTH+2)*MAX_IMAGE_BATCH*sizeof(float));
}

void CUDA_multiref_aligner::setup_params(vector<float> all_ali_params, vector<float> all_ctf_params) {

	ali_params = (float *)malloc(NIMA*4*sizeof(float));
	for (int i=0; i<NIMA*4; i++)   ali_params[i] = all_ali_params[i];
	if (CTF == 1) {
		ctf_params = (float *)malloc(NIMA*6*sizeof(float));
		for (int i=0; i<NIMA*6; i++)  ctf_params[i] = all_ctf_params[i];
	}
}

void CUDA_multiref_aligner::insert_image(EMData *image, int num) {

	int base_address = num*NX*NY;

	for (int y=0; y<NY; y++)
		for (int x=0; x<NX; x++)
			image_stack[base_address+y*NX+x] = (*image)(x, y);
}

void CUDA_multiref_aligner::insert_ref_image(EMData *image, int num) {

	int base_address = num*NX*NY;

	for (int y=0; y<NY; y++)
		for (int x=0; x<NX; x++)
			ref_image_stack[base_address+y*NX+x] = (*image)(x, y);
}

vector<float> CUDA_multiref_aligner::multiref_ali2d(int silent) {

	float *ctf_params_ref = (float *)malloc(NREF*6*sizeof(float));
	float *sx2 = (float *)malloc(NIMA*sizeof(float));
	float *sy2 = (float *)malloc(NIMA*sizeof(float));
	vector<float> align_results;
	int ccf_offset = NREF*(RING_LENGTH+2)*(2*KX+1)*(2*KY+1);

	vector<int> batch_size;
	vector<int> batch_begin;

	if (CTF == 1) {
		float previous_defocus = ctf_params[0];
		int current_size = 1;
		for (int i=1; i<NIMA; i++) {
			if (ctf_params[i*6] != previous_defocus || current_size >= MAX_IMAGE_BATCH) {
				batch_size.push_back(current_size);
				current_size = 1;
				previous_defocus = ctf_params[i*6];
			} else current_size++;
		}
		batch_size.push_back(current_size);
	} else {
		batch_size.resize(NIMA/MAX_IMAGE_BATCH, MAX_IMAGE_BATCH);
		if (NIMA%MAX_IMAGE_BATCH != 0)  batch_size.push_back(NIMA%MAX_IMAGE_BATCH);
	}
	int num_batch = batch_size.size();
	batch_begin.resize(num_batch, 0);
	for (int i=1; i<num_batch; i++) batch_begin[i] = batch_size[i-1]+batch_begin[i-1];
	assert(batch_begin[num_batch-1]+batch_size[num_batch-1] == NIMA-1);

	for (int i=0; i<NIMA; i++) {
		float ang = ali_params[i*4]/180.0*M_PI;
		float sx = ali_params[i*4+1];
		float sy = ali_params[i*4+2];
		float co = cos(ang);
		float so = sin(ang);
		sx2[i] = -(sx*co-sy*so);
		sy2[i] = -(sx*so+sy*co);
	}

	for (int i=0; i<num_batch; i++) {
		if (CTF == 1) {
			for (int p=0; p<NREF; p++)
				for (int q=0; q<6; q++)
					ctf_params_ref[p*6+q] = ctf_params[batch_begin[i]*6+q];
			filter_image(ref_image_stack, ref_image_stack_filtered, NREF, NX, NY, ctf_params_ref);
			calculate_multiref_ccf(image_stack+batch_begin[i]*NX*NY, ref_image_stack_filtered, ccf, batch_size[i], NREF, NX, NY, RING_LENGTH, NRING, OU, STEP, KX, KY,
				sx2+batch_begin[i], sy2+batch_begin[i], silent);
		} else {
			calculate_multiref_ccf(image_stack+batch_begin[i]*NX*NY, ref_image_stack, ccf, batch_size[i], NREF, NX, NY, RING_LENGTH, NRING, OU, STEP, KX, KY,
				sx2+batch_begin[i], sy2+batch_begin[i], silent);
		}

		for (int j=0; j<batch_size[i]; j++) {
			for (int im=0; im<NREF; im++) {
				float max_ccf = -1.0e22;
				float ang = 0.0, sx = 0.0, sy = 0.0;
				int mirror = 0;
				for (int kx=-KX; kx<=KX; kx++) {
					for (int ky=-KY; ky<=KY; ky++) {
						int base_address = (((ky+KY)*(2*KX+1)+(kx+KX))*NREF+im)*(RING_LENGTH+2)+ccf_offset*2*j;
						for (int l=0; l<RING_LENGTH; l++) {
							float ts = ccf[base_address+l];
							float tm = ccf[base_address+l+ccf_offset];
							if (ts > max_ccf) {
								ang = 360.0-float(l)/RING_LENGTH*360.0;
								sx = -kx*STEP;
								sy = -ky*STEP;
								mirror = 0;
								max_ccf = ts;
							}
							if (tm > max_ccf) {
								ang = float(l)/RING_LENGTH*360.0;
								sx = -kx*STEP;
								sy = -ky*STEP;
								mirror = 1;
								max_ccf = tm;
							}
						}
					}
				}
				float co =  cos(ang*M_PI/180.0);
				float so = -sin(ang*M_PI/180.0);

				int img_num = batch_begin[i]+j;
				float sxs = (sx-sx2[img_num])*co-(sy-sy2[img_num])*so;
				float sys = (sx-sx2[img_num])*so+(sy-sy2[img_num])*co;

				align_results.push_back(max_ccf);
				align_results.push_back(ang);
				align_results.push_back(sxs);
				align_results.push_back(sys);
				align_results.push_back(mirror);
			}
		}
	}

	free(ctf_params_ref);
	free(sx2);
	free(sy2);

	return align_results;
}

#endif


void EMAN::dump_aligners()
{
	dump_factory < Aligner > ();
}

map<string, vector<string> > EMAN::dump_aligners_list()
{
	return dump_factory_list < Aligner > ();
}
