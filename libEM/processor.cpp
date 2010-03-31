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
#include "plugins/processor_template.h"
#include "ctf.h"
#include "xydata.h"
#include "emdata.h"
#include "emassert.h"
#include "randnum.h"
#include "symmetry.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_wavelet2d.h>
#include <algorithm>
#include <ctime>

#ifdef EMAN2_USING_CUDA
#include "cuda/cuda_util.h"
#include "cuda/cuda_processor.h"
#endif // EMAN2_USING_CUDA

using namespace EMAN;
using std::reverse;

const string SNREvalProcessor::NAME = "eval.maskedsnr";
const string AmpweightFourierProcessor::NAME = "filter.ampweight";
const string ConvolutionProcessor::NAME = "math.convolution";
const string XGradientProcessor::NAME = "math.edge.xgradient";
const string YGradientProcessor::NAME = "math.edge.ygradient";
const string ZGradientProcessor::NAME = "math.edge.zgradient";
const string Wiener2DAutoAreaProcessor::NAME = "filter.wiener2dauto";
const string Wiener2DFourierProcessor::NAME = "filter.wiener2d";
const string LinearRampFourierProcessor::NAME = "filter.linearfourier";
const string LowpassSharpCutoffProcessor::NAME = "eman1.filter.lowpass.sharp";
const string HighpassSharpCutoffProcessor::NAME = "eman1.filter.highpass.sharp";
const string LowpassGaussProcessor::NAME = "eman1.filter.lowpass.gaussian";
const string HighpassAutoPeakProcessor::NAME = "filter.highpass.autopeak";
const string HighpassGaussProcessor::NAME = "eman1.filter.highpass.gaussian";
const string LowpassTanhProcessor::NAME = "eman1.filter.lowpass.tanh";
const string HighpassTanhProcessor::NAME = "eman1.filter.highpass.tanh";
const string HighpassButterworthProcessor::NAME = "eman1.filter.highpass.butterworth";
const string LinearRampProcessor::NAME = "eman1.filter.ramp";
const string AbsoluateValueProcessor::NAME = "math.absvalue";
const string BooleanProcessor::NAME = "threshold.notzero";
const string KmeansSegmentProcessor::NAME = "segment.kmeans";
const string InvertCarefullyProcessor::NAME = "math.invert.carefully";
const string ValuePowProcessor::NAME = "math.pow";
const string ValueSquaredProcessor::NAME = "math.squared";
const string ValueSqrtProcessor::NAME = "math.sqrt";
const string ToZeroProcessor::NAME = "threshold.belowtozero";
const string Rotate180Processor::NAME = "math.rotate.180";
const string TransformProcessor::NAME = "xform";
const string IntTranslateProcessor::NAME = "math.translate.int";
const string ScaleTransformProcessor::NAME = "xform.scale";
const string ClampingProcessor::NAME = "threshold.clampminmax";
const string NSigmaClampingProcessor::NAME = "threshold.clampminmax.nsigma";
const string ToMinvalProcessor::NAME = "threshold.belowtominval";
const string CutToZeroProcessor::NAME = "threshold.belowtozero_cut";
const string BinarizeProcessor::NAME = "threshold.binary";
const string BinarizeFourierProcessor::NAME = "threshold.binary.fourier";
const string CollapseProcessor::NAME = "threshold.compress";
const string LinearXformProcessor::NAME = "math.linear";
const string ExpProcessor::NAME = "math.exp";
const string FiniteProcessor::NAME = "math.finite";
const string RangeThresholdProcessor::NAME = "threshold.binaryrange";
const string SigmaProcessor::NAME = "math.sigma";
const string LogProcessor::NAME = "math.log";
const string MaskSharpProcessor::NAME = "mask.sharp";
const string MaskEdgeMeanProcessor::NAME = "mask.ringmean";
const string MaskNoiseProcessor::NAME = "mask.noise";
const string MaskGaussProcessor::NAME = "mask.gaussian";
const string MaskGaussNonuniformProcessor::NAME = "mask.gaussian.nonuniform";
const string MaskGaussInvProcessor::NAME = "math.gausskernelfix";
const string LinearPyramidProcessor::NAME = "math.linearpyramid";
const string MakeRadiusSquaredProcessor::NAME = "math.toradiussqr";
const string MakeRadiusProcessor::NAME = "math.toradius";
const string ComplexNormPixel::NAME = "complex.normpixels";
const string LaplacianProcessor::NAME = "math.laplacian";
const string ZeroConstantProcessor::NAME = "mask.contract";
const string BoxMedianProcessor::NAME = "eman1.filter.median";
const string BoxSigmaProcessor::NAME = "math.localsigma";
const string BoxMaxProcessor::NAME = "math.localmax";
const string MinusPeakProcessor::NAME = "math.submax";
const string PeakOnlyProcessor::NAME = "mask.onlypeaks";
const string DiffBlockProcessor::NAME = "eman1.filter.blockrange";
const string CutoffBlockProcessor::NAME = "eman1.filter.blockcutoff";
const string MaxShrinkProcessor::NAME = "math.maxshrink";
const string MinShrinkProcessor::NAME = "math.minshrink";
const string MeanShrinkProcessor::NAME = "math.meanshrink";
const string MedianShrinkProcessor::NAME = "math.medianshrink";
const string FFTResampleProcessor::NAME = "math.fft.resample";
const string GradientRemoverProcessor::NAME = "math.lineargradientfix";
const string GradientPlaneRemoverProcessor::NAME = "filter.gradientPlaneRemover";
const string FlattenBackgroundProcessor::NAME = "filter.flattenbackground";
const string RampProcessor::NAME = "filter.ramp";
const string VerticalStripeProcessor::NAME = "math.verticalstripefix";
const string RealToFFTProcessor::NAME = "math.realtofft";
const string SigmaZeroEdgeProcessor::NAME = "mask.zeroedgefill";
const string BeamstopProcessor::NAME = "mask.beamstop";
const string MeanZeroEdgeProcessor::NAME = "mask.dampedzeroedgefill";
const string AverageXProcessor::NAME = "math.averageovery";
const string DecayEdgeProcessor::NAME = "mask.decayedge2d";
const string ZeroEdgeRowProcessor::NAME = "mask.zeroedge2d";
const string ZeroEdgePlaneProcessor::NAME = "mask.zeroedge3d";
const string BilateralProcessor::NAME = "bilateral";
const string NormalizeUnitProcessor::NAME = "normalize.unitlen";
const string NormalizeUnitSumProcessor::NAME = "normalize.unitsum";
const string NormalizeStdProcessor::NAME = "normalize";
const string NormalizeMaskProcessor::NAME = "normalize.mask";
const string NormalizeRampNormVar::NAME = "normalize.ramp.normvar";
const string NormalizeByMassProcessor::NAME = "normalize.bymass";
const string NormalizeEdgeMeanProcessor::NAME = "normalize.edgemean";
const string NormalizeCircleMeanProcessor::NAME = "normalize.circlemean";
const string NormalizeLREdgeMeanProcessor::NAME = "normalize.lredge";
const string NormalizeMaxMinProcessor::NAME = "normalize.maxmin";
const string NormalizeRowProcessor::NAME = "normalize.rows";
const string NormalizeToLeastSquareProcessor::NAME = "normalize.toimage";
const string RotationalAverageProcessor::NAME = "math.rotationalaverage";
const string RotationalSubstractProcessor::NAME = "math.rotationalsubtract";
const string TransposeProcessor::NAME = "xform.transpose";
const string FlipProcessor::NAME = "xform.flip";
const string AddNoiseProcessor::NAME = "math.addnoise";
const string AddSigmaNoiseProcessor::NAME = "math.addsignoise";
const string AddRandomNoiseProcessor::NAME = "addspectralnoise";
const string FourierToCornerProcessor::NAME = "xform.fourierorigin.tocorner";
const string FourierToCenterProcessor::NAME = "xform.fourierorigin.tocenter";
const string PhaseToCenterProcessor::NAME = "xform.phaseorigin.tocenter";
const string PhaseToCornerProcessor::NAME = "xform.phaseorigin.tocorner";
const string AutoMask2DProcessor::NAME = "mask.auto2d";
const string AutoMaskAsymUnit::NAME = "mask.asymunit";
const string AutoMask3DProcessor::NAME = "mask.auto3d.thresh";
const string AutoMask3D2Processor::NAME = "mask.auto3d";
const string AddMaskShellProcessor::NAME = "mask.addshells";
const string PhaseToMassCenterProcessor::NAME = "xform.phasecenterofmass";
const string ToMassCenterProcessor::NAME = "xform.centerofmass";
const string ACFCenterProcessor::NAME = "xform.centeracf";
const string SNRProcessor::NAME = "eman1.filter.snr";
const string FileFourierProcessor::NAME = "eman1.filter.byfile";
const string SymSearchProcessor::NAME = "misc.symsearch";
const string LocalNormProcessor::NAME = "misc.localnorm";
const string IndexMaskFileProcessor::NAME = "mask.fromfile";
const string CoordinateMaskFileProcessor::NAME = "mask.fromfile.sizediff";
const string PaintProcessor::NAME = "mask.paint";
const string DirectionalSumProcessor::NAME = "misc.directional_sum";
const string WatershedProcessor::NAME = "watershed";
template<> const string BinaryOperateProcessor<MaxPixelOperator>::NAME = "operate.max";
template<> const string BinaryOperateProcessor<MinPixelOperator>::NAME = "operate.min";
const string MaxPixelOperator::NAME = "operate.max";
const string MinPixelOperator::NAME = "operate.min";
const string MatchSFProcessor::NAME = "filter.matchto";
const string SetSFProcessor::NAME = "filter.setstrucfac";
const string SmartMaskProcessor::NAME = "mask.smart";
const string IterBinMaskProcessor::NAME = "mask.addshells.gauss";
const string TestImagePureGaussian::NAME = "testimage.puregaussian";
const string TestImageFourierNoiseGaussian::NAME = "testimage.noise.fourier.gaussian";
const string TestImageFourierNoiseProfile::NAME = "testimage.noise.fourier.profile";
const string CTFSNRWeightProcessor::NAME = "ctf.snr.weight";
const string TestImageLineWave::NAME = "testimage.linewave";
const string TestTomoImage::NAME = "testimage.tomo.objects";
const string TestImageGradient::NAME = "testimage.gradient";
const string TestImageAxes::NAME = "testimage.axes";
const string TestImageGaussian::NAME = "testimage.gaussian";
const string TestImageScurve::NAME = "testimage.scurve";
const string TestImageSphericalWave::NAME = "testimage.sphericalwave";
const string TestImageSinewave::NAME = "testimage.sinewave";
const string TestImageSinewaveCircular::NAME = "testimage.sinewave.circular";
const string TestImageSquarecube::NAME = "testimage.squarecube";
const string TestImageEllipse::NAME = "testimage.ellipsoid";
const string TestImageHollowEllipse::NAME = "testimage.ellipsoid.hollow";
const string TestImageCirclesphere::NAME = "testimage.circlesphere";
const string TestImageNoiseUniformRand::NAME = "testimage.noise.uniform.rand";
const string TestImageNoiseGauss::NAME = "testimage.noise.gauss";
const string TestImageCylinder::NAME = "testimage.cylinder";
const string CCDNormProcessor::NAME = "filter.ccdnorm";
const string WaveletProcessor::NAME = "basis.wavelet";
const string TomoTiltEdgeMaskProcessor::NAME = "tomo.tiltedgemask";
const string TomoTiltAngleWeightProcessor::NAME = "tomo.tiltangleweight";
const string FFTProcessor::NAME = "basis.fft";
const string RadialProcessor::NAME = "mask.radialprofile";
const string HistogramBin::NAME = "histogram.bin";
const string ModelEMCylinderProcessor::NAME = "math.model_em_cylinder";
const string ApplyPolynomialProfileToHelix::NAME = "math.poly_radial_profile";
const string BinarySkeletonizerProcessor::NAME="gorgon.binary_skel";
const string MirrorProcessor::NAME = "mirror";
const string NewLowpassTopHatProcessor::NAME = "filter.lowpass.tophat";
const string NewHighpassTopHatProcessor::NAME = "filter.highpass.tophat";
const string NewBandpassTopHatProcessor::NAME = "filter.bandpass.tophat";
const string NewHomomorphicTopHatProcessor::NAME = "filter.homomorphic.tophat";
const string NewLowpassGaussProcessor::NAME = "filter.lowpass.gauss";
const string LowpassAutoBProcessor::NAME="filter.lowpass.autob";
const string NewHighpassGaussProcessor::NAME = "filter.highpass.gauss";
const string NewBandpassGaussProcessor::NAME = "filter.bandpass.gauss";
const string NewHomomorphicGaussProcessor::NAME = "filter.homomorphic.gauss";
const string NewInverseGaussProcessor::NAME = "filter.gaussinverse";
const string SHIFTProcessor::NAME = "filter.shift";
const string InverseKaiserI0Processor::NAME = "filter.kaiser_io_inverse";
const string InverseKaiserSinhProcessor::NAME = "filter.kaisersinhinverse";
const string NewRadialTableProcessor::NAME = "filter.radialtable";
const string NewLowpassButterworthProcessor::NAME = "filter.lowpass.butterworth";
const string NewHighpassButterworthProcessor::NAME = "filter.highpass.butterworth";
const string NewHomomorphicButterworthProcessor::NAME = "filter.homomorphic.butterworth";
const string NewLowpassTanhProcessor::NAME = "filter.lowpass.tanh";
const string NewHighpassTanhProcessor::NAME = "filter.highpass.tanh";
const string NewHomomorphicTanhProcessor::NAME = "filter.homomorphic.tanh";
const string NewBandpassTanhProcessor::NAME = "filter.bandpass.tanh";
const string CTF_Processor::NAME = "filter.CTF_";

#ifdef EMAN2_USING_CUDA
const string CudaMultProcessor::NAME = "cuda.math.mult";
const string CudaCorrelationProcessor::NAME = "cuda.correlate";
#endif //EMAN2_USING_CUDA

#if 0
//const string XYZProcessor::NAME = "XYZ";
#endif	//0


template <> Factory < Processor >::Factory()
{
	force_add< LowpassSharpCutoffProcessor >();
	force_add<HighpassSharpCutoffProcessor>();
	force_add<LowpassGaussProcessor>();
	force_add<HighpassGaussProcessor>();
	force_add<HighpassAutoPeakProcessor>();
	force_add<LinearRampFourierProcessor>();

	force_add<LowpassTanhProcessor>();
	force_add<HighpassTanhProcessor>();
	force_add<HighpassButterworthProcessor>();
	force_add<AmpweightFourierProcessor>();
	force_add<Wiener2DFourierProcessor>();
	force_add<LowpassAutoBProcessor>();

	force_add<LinearPyramidProcessor>();
	force_add<LinearRampProcessor>();
	force_add<AbsoluateValueProcessor>();
	force_add<BooleanProcessor>();
	force_add<KmeansSegmentProcessor>();
	force_add<ValuePowProcessor>();
	force_add<ValueSquaredProcessor>();
	force_add<ValueSqrtProcessor>();
	force_add<Rotate180Processor>();
	force_add<TransformProcessor>();
	force_add<ScaleTransformProcessor>();
	force_add<IntTranslateProcessor>();
	force_add<InvertCarefullyProcessor>();

	force_add<ClampingProcessor>();
	force_add<NSigmaClampingProcessor>();

	force_add<ToZeroProcessor>();
	force_add<ToMinvalProcessor>();
	force_add<CutToZeroProcessor>();
	force_add<BinarizeProcessor>();
	force_add<BinarizeFourierProcessor>();
	force_add<CollapseProcessor>();
	force_add<LinearXformProcessor>();

	force_add<ExpProcessor>();
	force_add<RangeThresholdProcessor>();
	force_add<SigmaProcessor>();
	force_add<LogProcessor>();
	force_add<FiniteProcessor>();

	force_add< BinaryOperateProcessor<MaxPixelOperator> >();
	force_add< BinaryOperateProcessor<MinPixelOperator> >();

	force_add<PaintProcessor>();
	force_add<WatershedProcessor>();
	force_add<MaskSharpProcessor>();
	force_add<MaskEdgeMeanProcessor>();
	force_add<MaskNoiseProcessor>();
	force_add<MaskGaussProcessor>();
	force_add<MaskGaussNonuniformProcessor>();
	force_add<MaskGaussInvProcessor>();

	force_add<MaxShrinkProcessor>();
	force_add<MinShrinkProcessor>();
	force_add<MeanShrinkProcessor>();
	force_add<MedianShrinkProcessor>();
	force_add<FFTResampleProcessor>();

	force_add<MakeRadiusSquaredProcessor>();
	force_add<MakeRadiusProcessor>();

	force_add<ComplexNormPixel>();

	force_add<LaplacianProcessor>();
	force_add<ZeroConstantProcessor>();

	force_add<BoxMedianProcessor>();
	force_add<BoxSigmaProcessor>();
	force_add<BoxMaxProcessor>();

	force_add<MinusPeakProcessor>();
	force_add<PeakOnlyProcessor>();
	force_add<DiffBlockProcessor>();

	force_add<CutoffBlockProcessor>();
	force_add<GradientRemoverProcessor>();
	force_add<GradientPlaneRemoverProcessor>();
	force_add<FlattenBackgroundProcessor>();
	force_add<VerticalStripeProcessor>();
	force_add<RealToFFTProcessor>();
	force_add<SigmaZeroEdgeProcessor>();
	force_add<RampProcessor>();

	force_add<BeamstopProcessor>();
	force_add<MeanZeroEdgeProcessor>();
	force_add<AverageXProcessor>();
	force_add<DecayEdgeProcessor>();
	force_add<ZeroEdgeRowProcessor>();
	force_add<ZeroEdgePlaneProcessor>();

	force_add<BilateralProcessor>();

	force_add<ConvolutionProcessor>();

	force_add<NormalizeStdProcessor>();
	force_add<NormalizeUnitProcessor>();
	force_add<NormalizeUnitSumProcessor>();
	force_add<NormalizeMaskProcessor>();
	force_add<NormalizeEdgeMeanProcessor>();
	force_add<NormalizeCircleMeanProcessor>();
	force_add<NormalizeLREdgeMeanProcessor>();
	force_add<NormalizeMaxMinProcessor>();
	force_add<NormalizeByMassProcessor>();
	force_add<NormalizeRowProcessor>();
	force_add<NormalizeRampNormVar>();

	force_add<HistogramBin>();

	force_add<NormalizeToLeastSquareProcessor>();

	force_add<RotationalAverageProcessor>();
	force_add<RotationalSubstractProcessor>();
	force_add<FlipProcessor>();
	force_add<TransposeProcessor>();
	force_add<MirrorProcessor>();

	force_add<AddNoiseProcessor>();
	force_add<AddSigmaNoiseProcessor>();
	force_add<AddRandomNoiseProcessor>();

	force_add<PhaseToCenterProcessor>();
	force_add<PhaseToCornerProcessor>();
	force_add<FourierToCenterProcessor>();
	force_add<FourierToCornerProcessor>();
	force_add<AutoMask2DProcessor>();
	force_add<AutoMask3DProcessor>();
	force_add<AutoMask3D2Processor>();
	force_add<AddMaskShellProcessor>();
	force_add<AutoMaskAsymUnit>();

	force_add<CTFSNRWeightProcessor>();

	force_add<ToMassCenterProcessor>();
	force_add<PhaseToMassCenterProcessor>();
	force_add<ACFCenterProcessor>();
	force_add<SNRProcessor>();

	force_add<XGradientProcessor>();
	force_add<YGradientProcessor>();
	force_add<ZGradientProcessor>();

	force_add<FileFourierProcessor>();

	force_add<SymSearchProcessor>();
	force_add<LocalNormProcessor>();

	force_add<IndexMaskFileProcessor>();
	force_add<CoordinateMaskFileProcessor>();
	force_add<SetSFProcessor>();
	force_add<MatchSFProcessor>();

	force_add<SmartMaskProcessor>();
	force_add<IterBinMaskProcessor>();

	force_add<TestImageGaussian>();
	force_add<TestImagePureGaussian>();
	force_add<TestImageSinewave>();
	force_add<TestImageSphericalWave>();
	force_add<TestImageSinewaveCircular>();
	force_add<TestImageSquarecube>();
	force_add<TestImageCirclesphere>();
	force_add<TestImageAxes>();
	force_add<TestImageNoiseUniformRand>();
	force_add<TestImageNoiseGauss>();
	force_add<TestImageScurve>();
	force_add<TestImageCylinder>();
	force_add<TestImageGradient>();
	force_add<TestTomoImage>();
	force_add<TestImageLineWave>();
	force_add<TestImageEllipse>();
	force_add<TestImageHollowEllipse>();
	force_add<TestImageFourierNoiseGaussian>();
	force_add<TestImageFourierNoiseProfile>();

	force_add<TomoTiltEdgeMaskProcessor>();
	force_add<TomoTiltAngleWeightProcessor>();

	force_add<NewLowpassTopHatProcessor>();
	force_add<NewHighpassTopHatProcessor>();
	force_add<NewBandpassTopHatProcessor>();
	force_add<NewHomomorphicTopHatProcessor>();
	force_add<NewLowpassGaussProcessor>();
	force_add<NewHighpassGaussProcessor>();
	force_add<NewBandpassGaussProcessor>();
	force_add<NewHomomorphicGaussProcessor>();
	force_add<NewInverseGaussProcessor>();
	force_add<NewLowpassButterworthProcessor>();
	force_add<NewHighpassButterworthProcessor>();
	force_add<NewHomomorphicButterworthProcessor>();
	force_add<NewLowpassTanhProcessor>();
	force_add<NewHighpassTanhProcessor>();
	force_add<NewBandpassTanhProcessor>();
	force_add<NewHomomorphicTanhProcessor>();
	force_add<NewRadialTableProcessor>();
	force_add<InverseKaiserI0Processor>();
	force_add<InverseKaiserSinhProcessor>();
	force_add<CCDNormProcessor>();
	force_add<CTF_Processor>();
	force_add<SHIFTProcessor>();

	force_add<WaveletProcessor>();
	force_add<FFTProcessor>();
	force_add<RadialProcessor>();

	force_add<DirectionalSumProcessor>();

	//Gorgon-related processors
	force_add<ModelEMCylinderProcessor>();
	force_add<ApplyPolynomialProfileToHelix>();
	force_add<BinarySkeletonizerProcessor>();

#ifdef EMAN2_USING_CUDA
	force_add<CudaMultProcessor>();
	force_add<CudaCorrelationProcessor>();
#endif // EMAN2_USING_CUDA
	
//	force_add<XYZProcessor>();
}

void FiniteProcessor::process_pixel(float *x) const
{
	if ( !Util::goodf(x) ) {
		*x = to;
	}
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

	preprocess(image);

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

		//ift->update(); Unecessary

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

void FourierAnlProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	preprocess(image);

// 	int array_size = FFTRADIALOVERSAMPLE * image->get_ysize();
// 	float step=0.5f/array_size;
//
// 	vector < float >yarray(array_size);


	if (image->is_complex()) {
		vector <float>yarray = image->calc_radial_dist(image->get_ysize()/2,0,1.0,1);
		create_radial_func(yarray,image);
		image->apply_radial_func(0, 0.5f/yarray.size(), yarray);
	}
	else {
		EMData *fft = image->do_fft();
		vector <float>yarray = fft->calc_radial_dist(fft->get_ysize()/2,0,1.0,1);
		create_radial_func(yarray,image);
		fft->apply_radial_func(0,  0.5f/yarray.size(), yarray);
		EMData *ift = fft->do_ift();

		memcpy(image->get_data(),ift->get_data(),ift->get_xsize()*ift->get_ysize()*ift->get_zsize()*sizeof(float));

		//ift->update(); Unecessary

		delete fft;
		delete ift;

	}

	image->update();
}

void LowpassAutoBProcessor::create_radial_func(vector < float >&radial_mask,EMData *image) const{
	float apix=(float)image->get_attr("apix_x");
	int verbose=(int)params["verbose"];
//	int adaptnoise=params.set_default("adaptnoise",0);
	float noisecutoff=(float)params.set_default("noisecutoff",0.0);
	if (apix<=0 || apix>7.0) throw ImageFormatException("0 < apix_x < 7.0");
	float ds=1.0/(apix*image->get_xsize());	// 0.5 is because radial mask is 2x oversampled 
	int start=(int)floor(1.0/(15.0*ds));
	int end=radial_mask.size()-2;
	if (noisecutoff>0) end=(int)floor(noisecutoff/ds);
	if (end>radial_mask.size()-2) {
		printf("WARNING: specified noisecutoff too close to Nyquist, reset !");
		end=radial_mask.size()-2;
	}
	if (end<start+2) {
		printf("WARNING: noise cutoff too close to 15 A ! Results will not be good...");
		start=end-5;
	}

	FILE *out=NULL;
	if (verbose>2)  out=fopen("fitplot.txt","w");
	int N=(radial_mask.size()-start-2);
	float *x=(float *)malloc(N*sizeof(float));
	float *y=(float *)malloc(N*sizeof(float));
	float *dy=(float *)malloc(N*sizeof(float));
	for (int i=start; i<radial_mask.size()-2; i++ ) {		// -2 is arbitrary because sometimes the last pixel or two have funny values
		x[i-start]=ds*ds*i*i;
		if (radial_mask[i]>0) y[i-start]=log(radial_mask[i]); // ok
		else if (i>start) y[i-start]=y[i-start-1];		// not good
		else y[i-start]=0.0;							// bad
		if (i<radial_mask.size()-3) dy[i-start]=y[i-start]-y[i-start-1];	// creates a 'derivative' of sorts, for use in adaptnoise
		if (out) fprintf(out,"%f\t%f\n",x[i-start],y[i-start]);
	}
	if (out) fclose(out);
	
	float slope=0,intercept=0;
	Util::calc_least_square_fit(end-start, x,y,&slope,&intercept,1);
 
	if (verbose) printf("slope=%f  intercept=%f\n",slope,intercept);
	
	float B=(float)params["bfactor"]+slope*4.0/2.0;		// *4 is for Henderson definition, 2.0 is for intensity vs amplitude
	float B2=(float)params["bfactor"];

	if (verbose) printf("User B = %1.2f  Corrective B = %1.2f  Total B=%1.3f\n",(float)params["bfactor"],slope*2.0,B);
	
	float cutval=exp(-B*pow(end*ds,2.0f)/4.0f)/exp(-B2*pow(end*ds,2.0f)/4.0f);
	for (int i=0; i<radial_mask.size(); i++) {
		if (i<=end) radial_mask[i]=exp(-B*pow(i*ds,2.0f)/4.0f);
		else radial_mask[i]=cutval*exp(-B2*pow(i*ds,2.0f)/4.0f);
	}
	if (verbose>1) Util::save_data(0,ds,radial_mask,"filter.txt");

	free(x);
	free(y);
	free(dy);
 }
	


void LowpassFourierProcessor::preprocess(EMData * image)
{
	if(params.has_key("apix")) {
		image->set_attr("apix_x", (float)params["apix"]);
		image->set_attr("apix_y", (float)params["apix"]);
		image->set_attr("apix_z", (float)params["apix"]);
	}

	const Dict dict = image->get_attr_dict();

	if( params.has_key("cutoff_abs") ) {
		lowpass = params["cutoff_abs"];
	}
	else if( params.has_key("cutoff_freq") ) {
		lowpass = (float)params["cutoff_freq"] * (float)dict["apix_x"] * (float)dict["nx"] / 2.0f;
	}
	else if( params.has_key("cutoff_pixels") ) {
		lowpass = (float)params["cutoff_pixels"] / (float)dict["nx"];
	}
}

void HighpassFourierProcessor::preprocess(EMData * image)
{
	if(params.has_key("apix")) {
		image->set_attr("apix_x", (float)params["apix"]);
		image->set_attr("apix_y", (float)params["apix"]);
		image->set_attr("apix_z", (float)params["apix"]);
	}

	const Dict dict = image->get_attr_dict();

	if( params.has_key("cutoff_abs") ) {
		highpass = params["cutoff_abs"];
	}
	else if( params.has_key("cutoff_freq") ) {
		highpass = (float)params["cutoff_freq"] * (float)dict["apix_x"] * (float)dict["nx"] / 2.0f;
	}
	else if( params.has_key("cutoff_pixels") ) {
		highpass = (float)params["cutoff_pixels"] / (float)dict["nx"];
	}
}

void SNREvalProcessor::process_inplace(EMData * image)
{
	int ys=image->get_ysize();

	EMData *mask1=new EMData(ys,ys,1);
	mask1->process_inplace("mask.gaussian",Dict("outer_radius", ys/2.0));
	EMData *mask2=mask1->copy();
	mask2->mult(-1.0f);
	mask2->add(1.0);
	mask2->process_inplace("mask.decayedge2d",Dict("width",4));

/*



mask1=EMData(ys2,ys2,1)
		mask1.to_one()
		mask1.process_inplace("mask.gaussian",{"outer_radius":radius})
		mask2=mask1.copy()*-1+1
#		mask1.process_inplace("mask.decayedge2d",{"width":4})
		mask2.process_inplace("mask.decayedge2d",{"width":4})
		mask1.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
		mask2.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
		ratio1=mask1.get_attr("square_sum")/(ys*ys)	#/1.035
		ratio2=mask2.get_attr("square_sum")/(ys*ys)
		masks[(ys,radius)]=(mask1,ratio1,mask2,ratio2)



*/
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
#ifdef	_WIN32
		else c = static_cast<float>(_hypot(fftd[i],fftd[i+1]));
#else
		else c = static_cast<float>(hypot(fftd[i],fftd[i+1]));
#endif	//_WIN32
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

EMData* KmeansSegmentProcessor::process(const EMData * const image)
{
	EMData * result = image->copy();
	
	int nseg = params.set_default("nseg",12);
	float thr = params.set_default("thr",-1.0e30f);
	int ampweight = params.set_default("ampweight",1);
	float maxsegsize = params.set_default("maxsegsize",10000.0f);
	float minsegsep = params.set_default("minsegsep",0.0f);
	int maxiter = params.set_default("maxiter",100);
	int maxvoxmove = params.set_default("maxvoxmove",25);
	int verbose = params.set_default("verbose",0);
	
	vector<float> centers(nseg*3);
	vector<float> count(nseg);
	int nx=image->get_xsize();
	int ny=image->get_ysize();
	int nz=image->get_zsize();
//	int nxy=nx*ny;
	
	// seed
	for (int i=0; i<nseg*3; i+=3) {
		centers[i]=  Util::get_frand(0.0f,(float)nx);
		centers[i+1]=Util::get_frand(0.0f,(float)ny);
		centers[i+2]=Util::get_frand(0.0f,(float)nz);
	}
	
	for (int iter=0; iter<maxiter; iter++) {
		// **** classify
		size_t pixmov=0;		// count of moved pixels
		for (int z=0; z<nz; z++) {
			for (int y=0; y<ny; y++) {
				for (int x=0; x<nz; x++) {
					if (image->get_value_at(x,y,z)<thr) {
						result->set_value_at(x,y,z,-1.0);		//below threshold -> -1 (unclassified)
						continue;
					}
					int bcls=-1;			// best matching class		
					float bdist=(float)(nx+ny+nz);	// distance for best class
					for (int c=0; c<nseg; c++) { 
						float d=Util::hypot3(x-centers[c*3],y-centers[c*3+1],z-centers[c*3+2]);
						if (d<bdist) { bdist=d; bcls=c; }
					}
					if ((int)result->get_value_at(x,y,z)!=bcls) pixmov++;
					if (bdist>maxsegsize) result->set_value_at(x,y,z,-1);		// pixel is too far from any center
					else result->set_value_at(x,y,z,(float)bcls);		// set the pixel to the class number
				}
			}
		}
	
		// **** adjust centers
		for (int i=0; i<nseg*3; i++) centers[i]=0;
		for (int i=0; i<nseg; i++) count[i]=0;
	
		// weighted sums
		for (int z=0; z<nz; z++) {
			for (int y=0; y<ny; y++) {
				for (int x=0; x<nz; x++) {
					int cls = (int)result->get_value_at(x,y,z);
					if (cls==-1) continue;
					float w=1.0;
					if (ampweight) w=image->get_value_at(x,y,z);
					
					centers[cls*3]+=x*w;
					centers[cls*3+1]+=y*w;
					centers[cls*3+2]+=z*w;
					count[cls]+=w;
				}
			}
		}
		
		// now each becomes center of mass, or gets randomly reseeded
		int nreseed=0;
		for (int c=0; c<nseg; c++) {
			// reseed
			if (count[c]==0) {
				nreseed++;
				do {
					centers[c*3]=  Util::get_frand(0.0f,(float)nx);
					centers[c*3+1]=Util::get_frand(0.0f,(float)ny);
					centers[c*3+2]=Util::get_frand(0.0f,(float)nz);
				} while (image->get_value_at((int)centers[c*3],(int)centers[c*3+1],(int)centers[c*3+2])<thr);		// This makes sure the new point is inside density
			}
			// center of mass
			else {
				centers[c*3]/=count[c];
				centers[c*3+1]/=count[c];
				centers[c*3+2]/=count[c];
			}
		}
		
		// with minsegsep, check separation
		if (minsegsep>0) {
			for (int c1=0; c1<nseg-1; c1++) {
				for (int c2=c1+1; c2<nseg; c2++) {
					if (Util::hypot3(centers[c1*3]-centers[c2*3],centers[c1*3+1]-centers[c2*3+1],centers[c1*3+2]-centers[c2*3+2])<=minsegsep) {
						nreseed++;
						do {
							centers[c1*3]=  Util::get_frand(0.0f,(float)nx);
							centers[c1*3+1]=Util::get_frand(0.0f,(float)ny);
							centers[c1*3+2]=Util::get_frand(0.0f,(float)nz);
						} while (image->get_value_at((int)centers[c1*3],(int)centers[c1*3+1],(int)centers[c1*3+2])<thr);
					}
				}
			}
		}

		
		if (verbose) printf("Iteration %3d: %6ld voxels moved, %3d classes reseeded\n",iter,pixmov,nreseed);
		if (nreseed==0 && pixmov<(size_t)maxvoxmove) break;		// termination conditions met
	}
	
	result->set_attr("segment_centers",centers);

	return result;
}

void KmeansSegmentProcessor::process_inplace(EMData *image)
{
	printf("Process inplace not implemented. Please use process.\n");
	return;
}


void LinearPyramidProcessor::process_inplace(EMData *image) {

	if (image->get_zsize()!=1) { throw ImageDimensionException("Only 2-D images supported"); }

	float *d=image->get_data();
	int nx=image->get_xsize();
	int ny=image->get_ysize();

	for (int y=0; y<ny; y++) {
		for (int x=0; x<nx; x++) {
			int l=x+y*nx;
			d[l]*=1.0f-abs(x-nx/2)*abs(y-ny/2)*4.0f/(nx*ny);
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


EMData * Wiener2DFourierProcessor::process(const EMData *in)
{
	const EMData *in2 = 0;
	if (in->is_complex()) in2=in;
	else in=in->do_fft();

	EMData *filt = in->copy_head();
	Ctf *ictf = ctf;

	if (!ictf) ctf=(Ctf *)in->get_attr("ctf");

	ictf->compute_2d_complex(filt,Ctf::CTF_WIENER_FILTER);
	filt->mult(*in2);
	EMData *ret=filt->do_ift();

	delete filt;
	if (!in->is_complex()) delete in2;

	if(!ictf && ctf) {delete ctf; ctf=0;}
	return(ret);
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

void LinearRampFourierProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	for (size_t i = 0; i < radial_mask.size(); i++) {
		radial_mask[i] = (float)i;
	}
}

void LowpassSharpCutoffProcessor::create_radial_func(vector < float >&radial_mask) const
{
	Assert(radial_mask.size() > 0);
	float x = 0.0f , step = 0.5f/radial_mask.size();
	printf("%d %f %f\n",(int)radial_mask.size(),lowpass,step);
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

void HighpassAutoPeakProcessor::preprocess(EMData * image)
{
	if(params.has_key("apix")) {
		image->set_attr("apix_x", (float)params["apix"]);
		image->set_attr("apix_y", (float)params["apix"]);
		image->set_attr("apix_z", (float)params["apix"]);
	}

	const Dict dict = image->get_attr_dict();

	if( params.has_key("cutoff_abs") ) {
		highpass = params["cutoff_abs"];
	}
	else if( params.has_key("cutoff_freq") ) {
		highpass = (float)params["cutoff_freq"] * (float)dict["apix_x"] * (float)dict["nx"] / 2.0f;
	}
	else if( params.has_key("cutoff_pixels") ) {
		highpass = (float)params["cutoff_pixels"] / (float)dict["nx"];
	}
}

void HighpassAutoPeakProcessor::create_radial_func(vector < float >&radial_mask, EMData *image) const
{
	unsigned int c;

//	for (unsigned int i=0; i<radial_mask.size(); i++) printf("%d\t%f\n",i,radial_mask[i]);
	for (c=2; c<radial_mask.size(); c++) if (radial_mask[c-1]<=radial_mask[c]) break;
	if (c>highpass) c=(unsigned int)highpass;		// the *2 is for the 2x oversampling

	radial_mask[0]=0.0;
//	for (int i=1; i<radial_mask.size(); i++) radial_mask[i]=(i<=c?radial_mask[c+1]/radial_mask[i]:1.0);
	for (unsigned int i=1; i<radial_mask.size(); i++) radial_mask[i]=(i<=c?0.0f:1.0f);

	printf("%f %d\n",highpass,c);
//	for (unsigned int i=0; i<radial_mask.size(); i++) printf("%d\t%f\n",i,radial_mask[i]);

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
	size_t i = 0;

	for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x++) {
				process_pixel(&data[i], x, y, z);
				++i;
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
		float r;
		for (int j=(y<r2?0:y-r2); j<(y+r2>ny?ny:y+r2); j++) {
			for (int i=(x<r2?0:x-r2); i<(x+r2>nx?nx:x+r2); i++) {
				r=sqrt(float(Util::square(i-x)+Util::square(j-y)));
				if (r>r2 && r>r1) continue;
				if (r>r1) image->set_value_at(i,j,0,v2*(r-r1)/(r2-r1)+v1*(r2-r)/(r2-r1));
				else image->set_value_at(i,j,0,v1);
			}
		}
	}
	else {
		float r;
		for (int k=(z<r2?0:z-r2); k<(z+r2>nz?nz:z+r2); k++) {
			for (int j=(y<r2?0:y-r2); j<(y+r2>ny?ny:y+r2); j++) {
				for (int i=(x<r2?0:x-r2); i<(x+r2>nx?nx:x+r2); i++) {
				r=sqrt(float(Util::square(i-x)+Util::square(j-y)+Util::square(k-z)));
				if (r>r2 && r>r1) continue;
				if (r>r1) image->set_value_at(i,j,k,v2*(r-r1)/(r2-r1)+v1*(r2-r)/(r2-r1));
				else image->set_value_at(i,j,k,v1);
				}
			}
		}
	}
	image->update();
}

void WatershedProcessor::process_inplace(EMData * image) {
	vector<float> xpoints = params["xpoints"];
	vector<float> ypoints = params["ypoints"];
	vector<float> zpoints = params["zpoints"];

	vector<int> x(xpoints.begin(),xpoints.end());
	vector<int> y(ypoints.begin(),ypoints.end());
	vector<int> z(zpoints.begin(),zpoints.end());


	// throw if vector lengths are unequal

//	float maxval = -99999;
	/*
	for(unsigned int i = 0; i < xpoints.size(); ++i) {
		float val = image->get_value_at(x[i],y[i],z[i]);
		if (val > maxval) {
			maxval = val;
		}
	}*/

	float minval = params["minval"];

	EMData* mask = new EMData(*image);
	mask->to_zero();

	// Set the original mask values
	for(unsigned int i = 0; i < xpoints.size(); ++i) {
		try {
			mask->set_value_at(x[i],y[i],z[i], (float)(i+1));
		} catch (...) {
			continue;
		}
	}
	mask->write_image("seeds2.mrc");
//	int dis = 500;
// 	float dx = (maxval-minval)/((float) dis - 1);


//	for(int i = 0; i < dis; ++i) {
//		float val = maxval-i*dx;

		while( true ) {
			bool cont= false;
			for(unsigned int j = 0; j < xpoints.size(); ++j)
			{

				Vec3i coord(x[j],y[j],z[j]);
				vector<Vec3i> region;
				region.push_back(coord);
				vector<Vec3i> find_region_input = region;
				while (true) {
					vector<Vec3i> v = find_region(mask,find_region_input, j+1, region);
					if (v.size() == 0 ) break;
					else find_region_input = v;
				}

				vector<Vec3i> tmp(region.begin(),region.end());
				region.clear();
				for(vector<Vec3i>::const_iterator it = tmp.begin(); it != tmp.end(); ++it ) {
					vector<Vec3i> tmp2 = watershed(mask, image, minval, *it, j+1);
					copy(tmp2.begin(),tmp2.end(),back_inserter(region));
				}
				if (region.size() != 0) cont = true;
			}

			if (!cont) break;
		}
//	}

	memcpy(image->get_data(),mask->get_data(),sizeof(float)*image->get_size());
	image->update();
}


vector<Vec3i > WatershedProcessor::find_region(EMData* mask,const vector<Vec3i >& coords, const int mask_value, vector<Vec3i >& region)
{
	static vector<Vec3i> two_six_connected;
	if (two_six_connected.size() == 0) {
		for(int i = -1; i <= 1; ++i) {
			for(int j = -1; j <= 1; ++j) {
				for(int  k = -1; k <= 1; ++k) {
					if ( j != 0 || i != 0 || k != 0) {
						two_six_connected.push_back(Vec3i(i,j,k));
					}
				}
			}
		}
	}

	vector<Vec3i> ret;
	for(vector<Vec3i>::const_iterator it = two_six_connected.begin(); it != two_six_connected.end(); ++it ) {
		for(vector<Vec3i>::const_iterator it2 = coords.begin(); it2 != coords.end(); ++it2 ) {
			if  (mask->get_value_at((*it2)[0],(*it2)[1],(*it2)[2]) != mask_value) throw;
			Vec3i c = (*it)+(*it2);

			if ( c[0] < 0 || c[0] >= mask->get_xsize()) continue;
			if ( c[1] < 0 || c[1] >= mask->get_ysize()) continue;
			if ( c[2] < 0 || c[2] >= mask->get_zsize()) continue;

			if( mask->get_value_at(c[0],c[1],c[2]) == mask_value ) {
				if (find(ret.begin(),ret.end(),c) == ret.end()) {
					if (find(region.begin(),region.end(),c) == region.end()) {
						region.push_back(c);
						ret.push_back(c);
					}
				}
			}
		}
	}
	return ret;
}

vector<Vec3i > WatershedProcessor::watershed(EMData* mask, EMData* image, const float& threshold, const Vec3i& coordinate, const int mask_value)
{
	static vector<Vec3i> two_six_connected;
	if (two_six_connected.size() == 0) {
		for(int i = -1; i <= 1; ++i) {
			for(int j = -1; j <= 1; ++j) {
				for(int  k = -1; k <= 1; ++k) {
					if ( j != 0 || i != 0 || k != 0) {
						two_six_connected.push_back(Vec3i(i,j,k));
					}
				}
			}
		}
	}

	if  (mask->get_value_at(coordinate[0],coordinate[1],coordinate[2]) != mask_value) throw;

	vector<Vec3i> ret;
	for(vector<Vec3i>::const_iterator it = two_six_connected.begin(); it != two_six_connected.end(); ++it ) {
		Vec3i c = (*it)+coordinate;

		if ( c[0] < 0 || c[0] >= image->get_xsize()) continue;
		if ( c[1] < 0 || c[1] >= image->get_ysize()) continue;
		if ( c[2] < 0 || c[2] >= image->get_zsize()) continue;

	//	cout << image->get_value_at(c[0],c[1],c[2] ) << " " << threshold << endl;
		if( image->get_value_at(c[0],c[1],c[2]) != 0 && (mask->get_value_at(c[0],c[1],c[2]) == 0 )) {
			//cout << "Added something " << mask_value << endl;
			mask->set_value_at(c[0],c[1],c[2], (float)mask_value);
			ret.push_back(c);
		}
	}
	return ret;
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
	size_t i = 0;

	for (int z = 0; z < nz; ++z) {
		for (int y = 0; y < ny; ++y) {
			for (int x = 0; x < nx; ++x) {
				float x1 = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc) + (z - zc) * (z - zc));
				if (x1 <= outer_radius + ring_width && x1 >= outer_radius - ring_width) {
					sum += data[i];
					++nitems;
				}
				++i;
			}
		}
	}

	ring_avg = sum / nitems;
}

void ToMinvalProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	float minval = params.set_default("minval",0.0f);
	float newval = params.set_default("newval",minval);

	size_t size = (size_t)image->get_xsize() *
		          (size_t)image->get_ysize() *
		          (size_t)image->get_zsize();
	float *data = image->get_data();

	

	for (size_t i = 0; i < size; i++) {
		if (data[i]<minval) data[i]=newval;
	}
	image->update();
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

	size_t k;
	for (int z = zstart; z < zend; z++) {
		for (int y = n; y < yend; y++) {
			for (int x = n; x < xend; x++) {

				k = z * nsec + y * nx + x;

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
		size_t knxy = k * nxy;

		for (int j = n; j < ny - n; j++) {
			int jnx = j * nx;

			for (int i = n; i < nx - n; i++) {
				size_t s = 0;

				for (int i2 = i - n; i2 <= i + n; i2++) {
					for (int j2 = j - n; j2 <= j + n; j2++) {
						for (int k2 = k - nzz; k2 <= k + nzz; k2++) {
							array[s] = data2[i2 + j2 * nx + k2 * nxy];
							++s;
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

#ifdef	_WIN32
					if (_hypot(j, i) < value2) {
#else
					if (hypot(j, i) < value2) {
#endif
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

void MedianShrinkProcessor::process_inplace(EMData * image)
{
	if (image->is_complex()) throw ImageFormatException("Error, the median shrink processor does not work on complex images");

	int shrink_factor =  params.set_default("n",0);
	if (shrink_factor <= 1) {
		throw InvalidValueException(shrink_factor,
									"median shrink: shrink factor must > 1");
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

// 	if ((nx % shrink_factor != 0) || (ny % shrink_factor != 0) || (nz > 1 && (nz % shrink_factor != 0))) {
// 		throw InvalidValueException(shrink_factor, "Image size not divisible by shrink factor");
// 	}


	int shrunken_nx = nx / shrink_factor;
	int shrunken_ny = ny / shrink_factor;
	int shrunken_nz = 1;
	if (nz > 1) shrunken_nz = nz / shrink_factor;

	EMData* copy = image->copy();
	image->set_size(shrunken_nx, shrunken_ny, shrunken_nz);
	accrue_median(image,copy,shrink_factor);
	image->update();
	if( copy )
	{
		delete copy;
		copy = 0;
	}
}

//
EMData* MedianShrinkProcessor::process(const EMData *const image)
{
	if (image->is_complex()) throw ImageFormatException("Error, the median shrink processor does not work on complex images");

	int shrink_factor =  params.set_default("n",0);
	if (shrink_factor <= 1) {
		throw InvalidValueException(shrink_factor,
									"median shrink: shrink factor must > 1");
	}
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();


// 	if ((nx % shrink_factor != 0) || (ny % shrink_factor != 0) || (nz > 1 && (nz % shrink_factor != 0))) {
// 		throw InvalidValueException(shrink_factor, "Image size not divisible by shrink factor");
// 	}


	int shrunken_nx = nx / shrink_factor;
	int shrunken_ny = ny / shrink_factor;
	int shrunken_nz = 1;
	if (nz > 1) shrunken_nz = nz / shrink_factor;

//	EMData* ret = new EMData(shrunken_nx, shrunken_ny, shrunken_nz);
	EMData *ret = image->copy_head();
	ret->set_size(shrunken_nx, shrunken_ny, shrunken_nz);

	accrue_median(ret,image,shrink_factor);
	ret->update();
	return ret;
}

void MedianShrinkProcessor::accrue_median(EMData* to, const EMData* const from,const int shrink_factor)
{

	int nx_old = from->get_xsize();
	int ny_old = from->get_ysize();

	int threed_shrink_factor = shrink_factor * shrink_factor;
	int z_shrink_factor = 1;
	if (from->get_zsize() > 1) {
		threed_shrink_factor *= shrink_factor;
		z_shrink_factor = shrink_factor;
	}

	float *mbuf = new float[threed_shrink_factor];


	int nxy_old = nx_old * ny_old;

	int nx = to->get_xsize();
	int ny = to->get_ysize();
	int nz = to->get_zsize();
	int nxy_new = nx * ny;

	float * rdata = to->get_data();
	const float *const data_copy = from->get_const_data();

	for (int l = 0; l < nz; l++) {
		int l_min = l * shrink_factor;
		int l_max = l * shrink_factor + z_shrink_factor;
		int cur_l = l * nxy_new;

		for (int j = 0; j < ny; j++) {
			int j_min = j * shrink_factor;
			int j_max = (j + 1) * shrink_factor;
			int cur_j = j * nx + cur_l;

			for (int i = 0; i < nx; i++) {
				int i_min = i * shrink_factor;
				int i_max = (i + 1) * shrink_factor;

				int k = 0;
				for (int l2 = l_min; l2 < l_max; l2++) {
					int cur_l2 = l2 * nxy_old;

					for (int j2 = j_min; j2 < j_max; j2++) {
						int cur_j2 = j2 * nx_old + cur_l2;

						for (int i2 = i_min; i2 < i_max; i2++) {
							mbuf[k] = data_copy[i2 + cur_j2];
							++k;
						}
					}
				}

				for (k = 0; k < threed_shrink_factor / 2 + 1; k++) {
					for (int i2 = k + 1; i2 < threed_shrink_factor; i2++) {
						if (mbuf[i2] < mbuf[k]) {
							float f = mbuf[i2];
							mbuf[i2] = mbuf[k];
							mbuf[k] = f;
						}
					}
				}

				rdata[i + cur_j] = mbuf[threed_shrink_factor / 2];
			}
		}
	}

	if( mbuf )
	{
		delete[]mbuf;
		mbuf = 0;
	}

	to->scale_pixel((float)shrink_factor);
}

EMData* FFTResampleProcessor::process(const EMData *const image)
{
	float sample_rate = params.set_default("n",0.0f);
	if (sample_rate <= 0.0F  )  {
		throw InvalidValueException(sample_rate,	"sample rate must be >0 ");
	}

	EMData* result;
	if (image->is_complex()) result = image->copy();
	else result = image->do_fft();
	fft_resample(result,image,sample_rate);
	// The image may have been padded - we should shift it so that the phase origin is where FFTW expects it
	result->update();
	result->scale_pixel(sample_rate);
	return result;
}

void FFTResampleProcessor::process_inplace(EMData * image)
{
	if (image->is_complex()) throw ImageFormatException("Error, the fft resampling processor does not work on complex images");


	float sample_rate = params.set_default("n",0.0f);
	if (sample_rate <= 0.0F  )  {
		throw InvalidValueException(sample_rate,	"sample rate (n) must be >0 ");
	}

	fft_resample(image,image,sample_rate);

	image->scale_pixel(sample_rate);
	image->update();


}

void FFTResampleProcessor::fft_resample(EMData* to, const EMData *const from, const float& sample_rate) {
	int nx = from->get_xsize();
	int ny = from->get_ysize();
	int nz = from->get_zsize();

	int new_nx = static_cast<int>( static_cast<float> (nx) / sample_rate);
	int new_ny = static_cast<int>( static_cast<float> (ny) / sample_rate);
	int new_nz = static_cast<int>( static_cast<float> (nz) / sample_rate);

	if (new_nx == 0) throw UnexpectedBehaviorException("The resample rate causes the pixel dimensions in the x direction to go to zero");
	if (new_ny == 0) new_ny = 1;
	if (new_nz == 0) new_nz = 1;

	int ndim = from->get_ndim();
	if ( ndim < 3 ) {
		new_nz = 1;
	}
	if ( ndim < 2 ) {
		new_ny = 1;
	}

	int fft_x_correction = 1;
	if (new_nx % 2 == 0) fft_x_correction = 2;

	int fft_y_correction = 0;
	if (ny != 1 && new_ny % 2 == 0 && ny % 2 == 1) fft_y_correction = 1;
	else if (ny != 1 && new_ny % 2 == 1 && ny % 2 == 0) fft_y_correction = -1;

	int fft_z_correction = 0;
	if (nz != 1 && new_nz % 2 == 0 && nz % 2 == 1) fft_z_correction = 1;
	else if (nz != 1 && new_nz % 2 == 1 && nz % 2 == 0) fft_z_correction = -1;

	if ( ! to->is_complex()) to->do_fft_inplace();

	if (ndim != 1) to->process_inplace("xform.fourierorigin.tocenter");

	Region clip(0,(ny-new_ny)/2-fft_y_correction,(nz-new_nz)/2-fft_z_correction,new_nx+fft_x_correction,new_ny,new_nz);
	to->clip_inplace(clip);

	if (fft_x_correction == 1) to->set_fftodd(true);
	else to->set_fftodd(false);

	if (ndim != 1) to->process_inplace("xform.fourierorigin.tocorner");

	to->do_ift_inplace();
	to->depad_corner();

}


EMData* MeanShrinkProcessor::process(const EMData *const image)
{
	if (image->is_complex()) throw ImageFormatException("Error, the mean shrink processor does not work on complex images");

	if (image->get_ndim() == 1) { throw ImageDimensionException("Error, mean shrink works only for 2D & 3D images"); }

	float shrink_factor0 = params.set_default("n",0.0f);
	int shrink_factor = int(shrink_factor0);
	if (shrink_factor0 <= 1.0F || ((shrink_factor0 != shrink_factor) && (shrink_factor0 != 1.5F) ) ) {
		throw InvalidValueException(shrink_factor0,
									"mean shrink: shrink factor must be >1 integer or 1.5");
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();


	// here handle the special averaging by 1.5 for 2D case
	if (shrink_factor0==1.5 ) {
		if (nz > 1 ) throw InvalidValueException(shrink_factor0, "mean shrink: only support 2D images for shrink factor = 1.5");

		int shrunken_nx = (int(nx / 1.5)+1)/2*2;	// make sure the output size is even
		int shrunken_ny = (int(ny / 1.5)+1)/2*2;
		EMData* result = new EMData(shrunken_nx,shrunken_ny,1);

		accrue_mean_one_p_five(result,image);
		result->update();

		return result;
	}

	int shrunken_nx = nx / shrink_factor;
	int shrunken_ny = ny / shrink_factor;
	int shrunken_nz = 1;

	if (nz > 1) {
		shrunken_nz = nz / shrink_factor;
	}

//	EMData* result = new EMData(shrunken_nx,shrunken_ny,shrunken_nz);
	EMData* result = image->copy_head();
	result->set_size(shrunken_nx,shrunken_ny,shrunken_nz);
	accrue_mean(result,image,shrink_factor);

	result->update();

	return result;
}

void MeanShrinkProcessor::process_inplace(EMData * image)
{
	if (image->is_complex()) throw ImageFormatException("Error, the mean shrink processor does not work on complex images");

	if (image->get_ndim() == 1) { throw ImageDimensionException("Error, mean shrink works only for 2D & 3D images"); }

	float shrink_factor0 = params.set_default("n",0.0f);
	int shrink_factor = int(shrink_factor0);
	if (shrink_factor0 <= 1.0F || ((shrink_factor0 != shrink_factor) && (shrink_factor0 != 1.5F) ) ) {
		throw InvalidValueException(shrink_factor0,
									"mean shrink: shrink factor must be >1 integer or 1.5");
	}

/*	if ((nx % shrink_factor != 0) || (ny % shrink_factor != 0) ||
	(nz > 1 && (nz % shrink_factor != 0))) {
	throw InvalidValueException(shrink_factor,
	"Image size not divisible by shrink factor");
}*/

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	// here handle the special averaging by 1.5 for 2D case
	if (shrink_factor0==1.5 ) {
		if (nz > 1 ) throw InvalidValueException(shrink_factor0, "mean shrink: only support 2D images for shrink factor = 1.5");

		int shrunken_nx = (int(nx / 1.5)+1)/2*2;	// make sure the output size is even
		int shrunken_ny = (int(ny / 1.5)+1)/2*2;

		EMData* orig = image->copy();
		image->set_size(shrunken_nx, shrunken_ny, 1);	// now nx = shrunken_nx, ny = shrunken_ny
		image->to_zero();

		accrue_mean_one_p_five(image,orig);

		if( orig ) {
			delete orig;
			orig = 0;
		}
		image->update();

		return;
	}

	accrue_mean(image,image,shrink_factor);

	int shrunken_nx = nx / shrink_factor;
	int shrunken_ny = ny / shrink_factor;
	int shrunken_nz = 1;
	if (nz > 1) shrunken_nz = nz / shrink_factor;

	image->update();
	image->set_size(shrunken_nx, shrunken_ny, shrunken_nz);
}

void MeanShrinkProcessor::accrue_mean(EMData* to, const EMData* const from,const int shrink_factor)
{
	const float * const data = from->get_const_data();
	float* rdata = to->get_data();

	int nx = from->get_xsize();
	int ny = from->get_ysize();
	int nz = from->get_zsize();
	int nxy = nx*ny;


	int shrunken_nx = nx / shrink_factor;
	int shrunken_ny = ny / shrink_factor;
	int shrunken_nz = 1;
	int shrunken_nxy = shrunken_nx * shrunken_ny;

	int normalize_shrink_factor = shrink_factor * shrink_factor;
	int z_shrink_factor = 1;

	if (nz > 1) {
		shrunken_nz = nz / shrink_factor;
		normalize_shrink_factor *= shrink_factor;
		z_shrink_factor = shrink_factor;
	}

	float invnormfactor = 1.0f/(float)normalize_shrink_factor;

	for (int k = 0; k < shrunken_nz; k++) {
		int k_min = k * shrink_factor;
		int k_max = k * shrink_factor + z_shrink_factor;
		size_t cur_k = k * shrunken_nxy;

		for (int j = 0; j < shrunken_ny; j++) {
			int j_min = j * shrink_factor;
			int j_max = j * shrink_factor + shrink_factor;
			size_t cur_j = j * shrunken_nx + cur_k;

			for (int i = 0; i < shrunken_nx; i++) {
				int i_min = i * shrink_factor;
				int i_max = i * shrink_factor + shrink_factor;

				float sum = 0;
				for (int kk = k_min; kk < k_max; kk++) {
					size_t cur_kk = kk * nxy;

					for (int jj = j_min; jj < j_max; jj++) {
						size_t cur_jj = jj * nx + cur_kk;
						for (int ii = i_min; ii < i_max; ii++) {
							sum += data[ii + cur_jj];
						}
					}
				}
				rdata[i + cur_j] = sum * invnormfactor;
			}
		}
	}
	to->scale_pixel((float)shrink_factor);
}


void MeanShrinkProcessor::accrue_mean_one_p_five(EMData* to, const EMData * const from)
{
	int nx0 = from->get_xsize(), ny0 = from->get_ysize();	// the original size

	int nx = to->get_xsize(), ny = to->get_ysize();

	float *data = to->get_data();
	const float * const data0 = from->get_const_data();

	for (int j = 0; j < ny; j++) {
		int jj = int(j * 1.5);
		float jw0 = 1.0F, jw1 = 0.5F;	// 3x3 -> 2x2, so each new pixel should have 2.25 of the old pixels
		if ( j%2 ) {
			jw0 = 0.5F;
			jw1 = 1.0F;
		}
		for (int i = 0; i < nx; i++) {
			int ii = int(i * 1.5);
			float iw0 = 1.0F, iw1 = 0.5F;
			float w = 0.0F;

			if ( i%2 ) {
				iw0 = 0.5F;
				iw1 = 1.0F;
			}
			if ( jj < ny0 ) {
				if ( ii < nx0 ) {
					data[j * nx + i] = data0[ jj * nx0 + ii ] * jw0 * iw0 ;
					w += jw0 * iw0 ;
					if ( ii+1 < nx0 ) {
						data[j * nx + i] += data0[ jj * nx0 + ii + 1] * jw0 * iw1;
						w += jw0 * iw1;
					}
				}
				if ( jj +1 < ny0 ) {
					if ( ii < nx0 ) {
						data[j * nx + i] += data0[ (jj+1) * nx0 + ii ] * jw1 * iw0;
						w += jw1 * iw0;
						if ( ii+1 < nx0 ) {
							data[j * nx + i] += data0[ (jj+1) * nx0 + ii + 1] * jw1 * iw1;
							w += jw1 * iw1;
						}
					}
				}
			}
			if ( w>0 ) data[j * nx + i] /= w;
		}
	}

	to->update();
	to->scale_pixel((float)1.5);
}

// This would have to be moved into the header if it were required in other source files
template<class LogicOp>
EMData* BooleanShrinkProcessor::process(const EMData *const image, Dict& params)
{
	// The basic idea of this code is to iterate through each pixel in the output image
	// determining its value by investigation a region of the input image

	if (!image) throw NullPointerException("Attempt to max shrink a null image");

	if (image->is_complex() ) throw ImageFormatException("Can not max shrink a complex image");


	int shrink = params.set_default("n",2);
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
void BooleanShrinkProcessor::process_inplace(EMData * image, Dict& params)
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

void FlattenBackgroundProcessor::process_inplace(EMData * image)
{

	EMData* mask = params.set_default("mask",(EMData*)0);
	int radius = params.set_default("radius",0);

	if (radius != 0 && mask != 0) throw InvalidParameterException("Error - the mask and radius parameters are mutually exclusive.");

	if (mask == 0 && radius == 0) throw InvalidParameterException("Error - you must specify either the mask or the radius parameter.");

	// If the radius isn't 0, then turn the mask into the thing we want...
	bool deletemask = false;
	if (radius != 0) {
		mask = new EMData;
		int n = image->get_ndim();
		if (n==1){
			mask->set_size(2*radius+1);
		} else if (n==2) {
			mask->set_size(2*radius+1,2*radius+1);
		}
		else /*n==3*/ {
			mask->set_size(2*radius+1,2*radius+1,2*radius+1);
		}
		// assuming default behavior is to make a circle/sphere with using the radius of the mask
		mask->process_inplace("testimage.circlesphere");
	}

	// Double check that that mask isn't too big
	int mnx = mask->get_xsize(); int mny = mask->get_ysize(); int mnz = mask->get_zsize();
	int nx = image->get_xsize(); int ny = image->get_ysize(); int nz = image->get_zsize();
	int nxc = nx+mnx; int nyc = ny+mny; int nzc = nz+mnz;
	if (nz == 1) nzc = 1; // Sanity check
	if (ny == 1) nyc = 1; // Sanity check

	if ( mnx > nx || mny > ny || mnz > nz)
		throw ImageDimensionException("Can not flatten using a mask that is larger than the image.");

	// Get the normalization factor
	float normfac = 0.0;
	for (int i=0; i<mask->get_xsize()*mask->get_ysize()*mask->get_zsize(); ++i){
		normfac += mask->get_value_at(i);
	}
	// If the sum is zero the user probably doesn't understand that determining a measure of the mean requires
	// strictly positive numbers. The user has specified a mask that consists entirely of zeros, or the mask
	// has a mean of zero.
	if (normfac == 0) throw InvalidParameterException("Error - the pixels in the mask sum to zero. This breaks the flattening procedure");
	normfac = 1.0f/normfac;

	// The mask can now be automatically resized to the dimensions of the image
//	bool undoclip = false;

	Region r;
	if (ny == 1) r = Region((mnx-nxc)/2,nxc);
	else if (nz == 1) r = Region((mnx-nxc)/2, (mny-nyc)/2,nxc,nyc);
	else r = Region((mnx-nxc)/2, (mny-nyc)/2,(mnz-nzc)/2,nxc,nyc,nzc);
	mask->clip_inplace(r,0);
//	undoclip = true;
//	if ( mnx < nx || mny < ny || mnz < nz) {
//		Region r((mnx-nx)/2, (mny-ny)/2,(mnz-nz)/2,nx,ny,nz);
//		mask->clip_inplace(r);
//		undoclip = true;
//	}

	Region r2;
	if (ny == 1) r2 = Region((nx-nxc)/2,nxc);
	else if (nz == 1) r2 = Region((nx-nxc)/2, (ny-nyc)/2,nxc,nyc);
	else r2 = Region((nx-nxc)/2, (ny-nyc)/2,(nz-nzc)/2,nxc,nyc,nzc);
	image->clip_inplace(r2,image->get_edge_mean());
	// Finally do the convolution
	EMData* m = image->convolute(mask);
	// Normalize so that m is truly the local mean
	m->mult(normfac);
	// Before we can subtract, the mean must be phase shifted
	m->process_inplace("xform.phaseorigin.tocenter");
	// Subtract the local mean
//	image->write_image("a.mrc");
//	m->write_image("b.mrc");
	image->sub(*m); // WE'RE DONE!
	delete m;

	if (deletemask) {
		delete mask;
	} else { // I clipped it inplace, so undo this clipping so the user gets back what the put in
		Region r;
		if (ny == 1) r = Region((nxc-mnx)/2,mnx);
		else if (nz == 1) r = Region((nxc-mnx)/2, (nyc-mny)/2,mnx,mny);
		else r = Region((nxc-mnx)/2, (nyc-mny)/2,(nzc-mnz)/2,mnx,mny,mnz);
		mask->clip_inplace(r);
	}

	Region r3;
	if (ny == 1) r3 = Region((nxc-nx)/2,nx);
	else if (nz == 1) r3 = Region((nxc-nx)/2, (nyc-ny)/2,nx,ny);
	else r3 = Region((nxc-nx)/2, (nyc-ny)/2,(nzc-nz)/2,nx,ny,nz);
	image->clip_inplace(r3);
//	if ( undoclip ) {
//		Region r((nx-mnx)/2, (ny-mny)/2, (nz-mnz)/2,mnx,mny,mnz);
//		mask->clip_inplace(r);
//	}

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

#ifdef	_WIN32
			int r = Util::round(_hypot((float) i - cenx, (float) j - ceny));
#else
			int r = Util::round(hypot((float) i - cenx, (float) j - ceny));
#endif	//_WIN32

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

	size_t idx;
	for (int z = 0; z < nz; z++) {
		for (int x = 0; x < nx; x++) {
			double sum = 0;
			for (int y = 0; y < ny; y++) {
				idx = x + y * nx + z * nxy;
				sum += data[idx];
			}
			float mean = (float) sum / ny;

			for (int y = 0; y < ny; y++) {
				idx = x + y * nx + z * nxy;
				data[idx] = mean;
			}
		}
	}

	image->update();
}

void DecayEdgeProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	if (image->get_zsize() > 1) throw ImageDimensionException("3D model not supported");

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	float *d = image->get_data();
	int width = params["width"];

	for (int i=0; i<width; i++) {
		float frac=i/(float)width;
		for (int j=0; j<nx; j++) {
			d[j+i*nx]*=frac;
			d[nx*ny-j-i*nx-1]*=frac;
		}
		for (int j=0; j<ny; j++) {
			d[j*nx+i]*=frac;
			d[nx*ny-j*nx-i-1]*=frac;
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

void NormalizeRampNormVar::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("cannot do normalization on NULL image");
		return;
	}

	if (image->is_complex()) {
		LOGWARN("cannot do normalization on complex image");
		return;
	}

	image->process_inplace( "filter.ramp" );
	int nx = image->get_xsize();
	EMData mask(nx,nx);
	mask.process_inplace("testimage.circlesphere", Dict("radius",nx/2-2,"fill",1));

	vector<float> rstls = Util::infomask( image, &mask, false);
	image->add((float)-rstls[0]);
	image->mult((float)1.0/rstls[1]);
	image->update();
}

void NormalizeByMassProcessor::process_inplace(EMData * image)
{
	float mass = params.set_default("mass",-1.0f);

	if (mass <= 0) throw InvalidParameterException("You must specify a positive non zero mass");

	float thr = params.set_default("thr",(float)image->get_attr("mean")+(float)image->get_attr("sigma"));

	float apix = params.set_default("apix",-1.123456789f);
	if (apix == -1.123456789 ) {
		if (image->has_attr("apix_x")) {
			apix = image->get_attr("apix_x");
		}
	}

	if (apix <= 0) throw InvalidParameterException("You must specify a positive non zero apix");

	float step = ((float)image->get_attr("sigma"))/2.0f;

	int count=0;
	size_t n = image->get_size();
	float* d = image->get_data();

	for (size_t i=0; i<n; ++i) {
		if (d[i]>=thr) ++count;
	}

	float max = image->get_attr("maximum");
	float min = image->get_attr("minimum");
	for (int j=0; j<4; j++) {
		while (thr<max && count*apix*apix*apix*.81/1000.0>mass) {
			thr+=step;
			count=0;
			for (size_t i=0; i<n; ++i) {
				if (d[i]>=thr) ++count;
			}
		}

		step/=4.0;

		while (thr>min && count*apix*apix*apix*.81/1000.0<mass) {
			thr-=step;
			count=0;
			for (size_t i=0; i<n; ++i) {
				if (d[i]>=thr) ++count;
			}
		}

		step/=4.0;
	}

	image->mult((float)1.0/thr);
	image->update();
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
		size_t l = i * nx;
		size_t r = l + nx - 2;
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
	size_t size = nx * ny * nz;

	float sum_x = 0;
	float sum_y = 0;
	int count = 0;

	for (size_t i = 0; i < size; i++) {
		if (refp[i] >= low_threshold && refp[i] <= high_threshold && refp[i] != 0.0 && rawp[i] != 0.0) {
			count++;
			sum_x += refp[i];
			sum_y += rawp[i];
		}
	}

	float sum_x_mean = sum_x / count;
	float sum_tt = 0;
	float b = 0;

	float t;
	for (size_t i = 0; i < size; i++) {
		if (refp[i] >= low_threshold && refp[i] <= high_threshold && refp[i] != 0.0 && rawp[i] != 0.0) {
			t = refp[i] - sum_x_mean;
			sum_tt += t * t;
			b += t * rawp[i];
		}
	}

	b /= sum_tt;

	float a = (sum_y - sum_x * b) / count;
	float scale = 1 / b;
	float shift = -a / b;

	for (size_t i = 0; i < size; i++) {
		rawp[i] = (rawp[i] - a) / b;
	}

	image->update();

	params["scale"] = scale;
	params["shift"] = shift;

	image->set_attr("norm_mult",scale);
	image->set_attr("norm_add",shift);

}


void BinarizeFourierProcessor::process_inplace(EMData* image) {
	ENTERFUNC;
	if (!image->is_complex()) throw ImageFormatException("Fourier binary thresholding processor only works for complex images");

	float threshold = params.set_default("value",-1.0f);
	if (threshold < 0) throw InvalidParameterException("For fourier amplitude-based thresholding, the threshold must be greater than or equal to 0.");

	image->ri2ap(); //  works for cuda

#ifdef EMAN2_USING_CUDA
	if (image->gpu_operation_preferred()) {
		EMDataForCuda tmp = image->get_data_struct_for_cuda();
		binarize_fourier_amp_processor(&tmp,threshold);
		image->set_ri(true); // So it can be used for fourier multiplaction, for example
		image->gpu_update();
		EXITFUNC;
		return;
	}
#endif

	float* d = image->get_data();
	for( size_t i = 0; i < image->get_size()/2; ++i, d+=2) {
		float v = *d;
		if ( v >= threshold ) {
			*d = 1;
			*(d+1) = 0;
		} else {
			*d = 0;
			*(d+1) = 0;
		}
	}

	// No need to run ap2ri, because 1+0i is the same in either notation
	image->set_ri(true); // So it can be used for fourier multiplaction, for example
	image->update();
	EXITFUNC;
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

			size_t idx;
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

								idx = k * height * width + i * width + j;
								new_img[idx] = f2 / f1;
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

void RotationalAverageProcessor::process_inplace(EMData * image)
{
	if (!image || image->is_complex()) {
		LOGWARN("only works on real image. do nothing.");
		return;
	}

	if (image->get_ndim() <= 0 || image->get_ndim() > 3)	throw ImageDimensionException("radial average processor only works for 2D and 3D images");

	float *rdata = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();

	vector < float >dist = image->calc_radial_dist(nx / 2, 0, 1,0);

	float midx = (float)((int)nx/2);
	float midy = (float)((int)ny/2);

	size_t c = 0;
	if (image->get_ndim() == 2) {
		for (int y = 0; y < ny; y++) {
			for (int x = 0; x < nx; x++, c++) {
	#ifdef	_WIN32
				float r = (float) _hypot(x - midx, y - midy);
	#else
				float r = (float) hypot(x - midx, y - midy);
	#endif	//_WIN32


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
	}
	else if (image->get_ndim() == 3) {
		int nz = image->get_zsize();
		float midz = (float)((int)nz/2);
		float r;
		int i;
		for (int z = 0; z < nz; ++z) {
			for (int y = 0; y < ny; ++y) {
				for (int x = 0; x < nx; ++x, ++c) {

					r = (float) Util::hypot3(x - midx, y - midy, z - midz);

					i = Util::fast_floor(r);
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
		}
	}

	image->update();
}



void RotationalSubstractProcessor::process_inplace(EMData * image)
{
	if (!image || image->is_complex()) {
		LOGWARN("only works on real image. do nothing.");
		return;
	}

	if (image->get_ndim() != 2) throw ImageDimensionException("This processor works only for 2D images");

	float *rdata = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();

	vector < float >dist = image->calc_radial_dist(nx / 2, 0, 1,0);

	int c = 0;
	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++, c++) {
#ifdef	_WIN32
			float r = (float) _hypot(x - nx / 2, y - ny / 2);
#else
			float r = (float) hypot(x - nx / 2, y - ny / 2);
#endif
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


EMData* TransposeProcessor::process(const EMData* const image) {
	if (image->get_ndim() != 2) throw UnexpectedBehaviorException("Transpose processor only works with 2D images");
	if (image->is_complex()) throw UnexpectedBehaviorException("Transpose processor only works with real images");

	EMData* ret = new EMData(image->get_ysize(),image->get_xsize(),1); // transpose dimensions

	for(int j = 0; j< image->get_ysize();++j) {
		for(int i = 0; i< image->get_xsize();++i) {
			ret->set_value_at(j,i,image->get_value_at(i,j));
		}
	}

	return ret;

}

void TransposeProcessor::process_inplace(EMData* image) {
	if (image->get_ndim() != 2) throw UnexpectedBehaviorException("Transpose processor only works with 2D images");
	if (image->is_complex()) throw UnexpectedBehaviorException("Transpose processor only works with real images");

	float* data = (float*)malloc(image->get_ysize()*image->get_xsize()*sizeof(float));

	int nx = image->get_ysize(); // note tranpose
	for(int j = 0; j< image->get_ysize();++j) {
		for(int i = 0; i< image->get_xsize();++i) {
			data[i*nx+j] = image->get_value_at(i,j);
		}
	}

	image->set_data(data,image->get_ysize(),image->get_xsize(),1);

}

void FlipProcessor::process_inplace(EMData * image)
{
	ENTERFUNC;
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	string axis = (const char*)params["axis"];
#ifdef EMAN2_USING_CUDA
	if (image->gpu_operation_preferred()) {
		float array[12] = {1.0, 0.0, 0.0, 0.0,
						   0.0, 1.0, 0.0, 0.0,
		 				   0.0, 0.0, 1.0, 0.0};
		if (axis == "x" || axis == "X") {		// horizontal flip
			array[0] = -1.0;
		}else if (axis == "y" || axis == "Y") {		// vertical flip
			array[5] = -1.0;
		}
		else if (axis == "z" || axis == "Z") {		// vertical flip
			array[10] = -1.0;
		}
		Transform t(array);
		Dict params("transform",(Transform*)&t);
		image->process_inplace("xform",params);
		EXITFUNC;
		return;
	}
#endif


	float *d = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	size_t nxy = nx * ny;


	// Note in all cases the origin is nx/2, ny/2 and nz/2
	// This means when flipping even sized dimensions that some pixels are redundant.
	// Here redundant pixels are set to zero, however, should this change to something
	// like the mean.
	if (axis == "x" || axis == "X") {		// Horizontal flip
		int offset = (nx%2 == 0);
		size_t idx1, idx2;
		for(int z = 0; z < nz; ++z) {
			for(int y = 0; y < ny; ++y) {
				if (offset != 0 ) {
					idx1 = z*nxy + y*nx;
					d[idx1] = 0; // Here's where you'd make it the mean
				}
				for(int x = offset; x < nx / 2; ++x) {
					idx1 = z*nxy + y*nx + x;
					idx2 = z*nxy + y*nx + (nx-x-1+offset);
					std::swap(d[idx1], d[idx2]);
				}

			}
		}
	}

	else if (axis == "y" || axis == "Y") {		// vertical flip
		int offset = (ny%2 == 0);
		for(int z=0; z<nz; ++z) {
			if (offset != 0) {
				std::fill(d+z*nxy,d+z*nxy+nx,0); // So if we have change it to the mean it's easy to do so. (instead of using memset)
			}
			for(int y=offset; y<ny/2; ++y) {
				for(int x=0; x<nx; ++x) {
					std::swap(d[z*nxy + y*nx +x], d[z*nxy + (ny -y -1+offset)*nx +x]);
				}
			}
		}
	}
	else if (axis == "z" || axis == "Z") {		//z axis flip
		int offset = (nz%2 == 0);
		if (offset != 0) {
			std::fill(d,d+nxy,0);// So if we have change it to the mean it's easy to do so. (instead of using memset)
		}
		size_t idx1, idx2;
		for(int z=offset; z<nz/2; ++z) {
			for(int y=0; y<ny; ++y) {
				for(int x=0; x<nx; ++x) {
					idx1 = z*nxy + y*nx + x;
					idx2 = (nz-z-1+offset)*nxy + y*nx + x;
					std::swap(d[idx1], d[idx2]);
				}
			}
		}
	}

	image->update();
	EXITFUNC;
}

void AddNoiseProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	Randnum * randnum = Randnum::Instance();
	if(params.has_key("seed")) {
		randnum->set_seed((int)params["seed"]);
	}

	float addnoise = params["noise"];
	addnoise *= get_sigma(image);
	float *dat = image->get_data();

	for (size_t j = 0; j < image->get_size(); ++j) {
		dat[j] += randnum->get_gauss_rand(addnoise, addnoise / 2);
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
		size_t idx;
		for( int s = 0; s < nz; s++ ) {
			for( int c =0; c < nx; c += 2 ) {
				idx = s*nxy+ny/2*nx+c;
				prev[0] = rdata[idx];
				prev[1] = rdata[idx+1];
				for( int r = 0; r <= ny/2; ++r ) {
					idx = s*nxy+r*nx+c;
					float* p1 = &rdata[idx];
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
	size_t idx1, idx2;
	for( int s = 0; s < nz; ++s ) {
		for( int r = 0 + yodd; r < ny/2+yodd; ++r ) {
			for( int c =0; c < nx; c += 2 ) {
				idx1 = s*nxy+r*nx+c;
				idx2 = s*nxy+(r+ny/2)*nx+c;
				p1 = &rdata[idx1];
				p2 = &rdata[idx2];

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
			size_t idx;
			for( int r = 0; r < ny; ++r ) {
				for( int c =0; c < nx; c += 2 ) {
					idx = nz/2*nxy+r*nx+c;
					prev[0] = rdata[idx];
					prev[1] = rdata[idx+1];
					for( int s = 0; s <= nz/2; ++s ) {
						idx = s*nxy+r*nx+c;
						float* p1 = &rdata[idx];
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
		size_t idx1, idx2;
		for( int s = 0+zodd; s < nz/2 + zodd; ++s ) {
			for( int r = 0; r < ny; ++r ) {
				for( int c =0; c < nx; c += 2 ) {
					idx1 = s*nxy+r*nx+c;
					idx2 = (s+nz/2)*nxy+r*nx+c;
					p1 = &rdata[idx1];
					p2 = &rdata[idx2];

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
	image->set_shuffled(false);
}

void FourierToCenterProcessor::process_inplace(EMData * image)
{
//	if ( !image->is_complex() ) throw ImageFormatException("Can not Fourier origin shift an image that is not complex");

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

	// This will tackle the 'normalization' images which come out of the Fourier reconstructor. 
	// ie- real-space 1/2 FFt images centered on the corner
	if ( !image->is_complex() ) {
		if (nz!=1 && !yodd && !zodd) {
			for (int x=0; x<nx; x++) {
				for (int y=0; y<ny; y++) {
					for (int z=0; z<nz/2; z++) {
						int y2=(y+ny/2)%ny;
						int z2=(z+nz/2)%nz;		// %nz should be redundant here
						size_t i=x+y*nx+z*nxy;
						size_t i2=x+y2*nx+z2*nxy;
						float swp=rdata[i];
						rdata[i]=rdata[i2];
						rdata[i2]=swp;
					}
				}
			}
			
			return;
		}
		else throw ImageFormatException("Can not Fourier origin shift an image that is not complex unless it is even in ny,nz and nx=ny/2+1");
	}

	if (yodd){
		// In 3D this is swapping the bottom slice (with respect to the y direction) and the middle slice,
		// shifting all slices below the middle slice down one. In 2D it is equivalent, but in terms of rows.
		float prev[2];
		size_t idx;
		for( int s = 0; s < nz; s++ ) {
			for( int c =0; c < nx; c += 2 ) {
				idx = s*nxy+c;
				prev[0] = rdata[idx];
				prev[1] = rdata[idx+1];
				for( int r = ny/2; r >= 0; --r ) {
					idx = s*nxy+r*nx+c;
					float* p1 = &rdata[idx];
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
	size_t idx1, idx2;
	for( int s = 0; s < nz; ++s ) {
		for( int r = 0; r < ny/2; ++r ) {
			for( int c =0; c < nx; c += 2 ) {
				idx1 = s*nxy+r*nx+c;
				idx2 = s*nxy+(r+ny/2+yodd)*nx+c;
				p1 = &rdata[idx1];
				p2 = &rdata[idx2];

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
			size_t idx;
			for( int r = 0; r < ny; ++r ) {
				for( int c =0; c < nx; c += 2 ) {
					prev[0] = rdata[r*nx+c];
					prev[1] = rdata[r*nx+c+1];
					for( int s = nz/2; s >= 0; --s ) {
						idx = s*nxy+r*nx+c;
						float* p1 = &rdata[idx];
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
		size_t idx1, idx2;
		for( int s = 0; s < nz/2; ++s ) {
			for( int r = 0; r < ny; ++r ) {
				for( int c =0; c < nx; c += 2 ) {
					idx1 = s*nxy+r*nx+c;
					idx2 = (s+nz/2+zodd)*nxy+r*nx+c;
					p1 = &rdata[idx1];
					p2 = &rdata[idx2];

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
	image->set_shuffled(true);
}

void Phase180Processor::fourier_phaseshift180(EMData * image)
{
	if ( !image->is_complex() ) throw ImageFormatException("Can not handle images that are not complex in fourier phase shift 180");

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	int nxy = nx * ny;

	float *rdata = image->get_data();

	// Who uses this function? It doesn't work for odd images, and it will give incorrect results for some even images
	// d.woolford, March 15 2009
	int of=0;
	if (((ny/2)%2)+((nz/2)%2)==1) of=1;

	for (int k = 0; k < nz; k++) {
		size_t k2 = k * nxy;

		for (int j = 0; j < ny; j++) {
			int i = ((k+j)%2==of?2:0);
			size_t j2 = j * nx + k2;

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
		size_t idx1, idx2;

		for ( int s = 0; s < nz/2; ++s ) {
			for ( int r = 0; r < ny/2; ++r ) {
				for ( int c = 0; c < nx/2; ++ c) {
					idx1 = s*nxy+r*nx+c;
					idx2 = (s+nz/2+zodd)*nxy+(r+ny/2+yodd)*nx+c+nx/2+xodd;
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
					idx1 = s*nxy+r*nx+c;
					idx2 = (s+nz/2+zodd)*nxy+(r+ny/2+yodd)*nx+c-nx/2-xodd;
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
					idx1 = s*nxy+r*nx+c;
					idx2 = (s+nz/2+zodd)*nxy+(r-ny/2-yodd)*nx+c-nx/2-xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}
		}
		// Swap the top left front and back right bottom quadrants
		for ( int s = 0; s < nz/2; ++s ) {
			for ( int r = ny-1; r >= (ny/2+yodd); --r ) {
				for ( int c = 0; c < nx/2; ++c) {
					idx1 = s*nxy+r*nx+c;
					idx2 = (s+nz/2+zodd)*nxy+(r-ny/2-yodd)*nx+c+nx/2+xodd;
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
			size_t idx1, idx2;
			for( int s = 0; s < nz/2; ++s ) {
				for ( int r = 0; r < ny/2; ++r ) {
					idx1 = s*nxy+r*nx+c;
					idx2 = (s+nz/2+zodd)*nxy+(r+ny/2+yodd)*nx+c;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}

			for( int s = nz-1; s >= (nz/2+zodd); --s ) {
				for ( int r = 0; r < ny/2; ++r ) {
					idx1 = s*nxy+r*nx+c;
					idx2 = (s-nz/2-zodd)*nxy+(r+ny/2+yodd)*nx+c;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}
		}
		if ( yodd ) {
			// Iterate along the y = ny/2 slice, swapping values where appropriate
			int r = ny/2;
			size_t idx1, idx2;
			for( int s = 0; s < nz/2; ++s ) {
				for ( int c = 0; c < nx/2; ++c ) {
					idx1 = s*nxy+r*nx+c;
					idx2 =(s+nz/2+zodd)*nxy+r*nx+c+nx/2+xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}

			for( int s = nz-1; s >= (nz/2+zodd); --s ) {
				for ( int c = 0; c < nx/2; ++c ) {
					idx1 = s*nxy+r*nx+c;
					idx2 = (s-nz/2-zodd)*nxy+r*nx+c+nx/2+xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}
		}
		if ( zodd ) {
			// Iterate along the z = nz/2 slice, swapping values where appropriate
			int s = nz/2;
			size_t idx1, idx2;
			for( int r = 0; r < ny/2; ++r ) {
				for ( int c = 0; c < nx/2; ++c ) {
					idx1 = s*nxy+r*nx+c;
					idx2 = s*nxy+(r+ny/2+yodd)*nx+c+nx/2+xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}

			for( int r = ny-1; r >= (ny/2+yodd); --r ) {
				for ( int c = 0; c < nx/2; ++c ) {
					idx1 = s*nxy+r*nx+c;
					idx2 = s*nxy+(r-ny/2-yodd)*nx+c+nx/2+xodd;
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
			size_t idx = 0;
			for (int r = 0; r < ny; ++r){
				for (int c = 0; c < nx; ++c) {
					float last_val = rdata[(nz-1)*nxy+r*nx+c];
					for (int s = nz/2; s < nz; ++s) {
						idx = s*nxy+r*nx+c;
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
			size_t idx = 0;
			for (int s = 0; s < nz; ++s) {
				for (int c = 0; c < nx; ++c) {
				float last_val = rdata[s*nxy+(ny-1)*nx+c];
					for (int r = ny/2; r < ny; ++r){
						idx = s*nxy+r*nx+c;
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
			size_t idx = 0;
			for (int s = 0; s < nz; ++s) {
				for (int r = 0; r < ny; ++r) {
					float last_val = rdata[s*nxy+r*nx+nx-1];
					for (int c = nx/2; c < nx; ++c){
						idx = s*nxy+r*nx+c;
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
	bool proceed = true;

#ifdef EMAN2_USING_CUDA
	bool cpu = image->cpu_rw_is_current();
	bool gpu = image->gpu_rw_is_current();
	if ( !cpu && !gpu )
		throw UnexpectedBehaviorException("Both the CPU and GPU data are not current");
	if (gpu && image->get_ndim() == 2) { // Because CUDA phase origin to center only works for 2D atm
		EMDataForCuda tmp = image->get_data_struct_for_cuda();
		emdata_phaseorigin_to_center(&tmp);
		proceed = false;
		image->gpu_update();
	}
#endif // EMAN2_USING_CUDA
	if (!proceed) return; // GPU processing occurred

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
			size_t idx = 0;
			for (int s = 0; s < nz; ++s) {
				for (int r = 0; r < ny; ++r) {
					float last_val = rdata[s*nxy+r*nx+nx/2];
					for (int c = nx-1; c >= nx/2; --c){
						idx = s*nxy+r*nx+c;
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
			size_t idx = 0;
			for (int s = 0; s < nz; ++s) {
				for (int c = 0; c < nx; ++c) {
					float last_val = rdata[s*nxy+ny/2*nx+c];
					for (int r = ny-1; r >= ny/2; --r){
						idx = s*nxy+r*nx+c;
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
			size_t idx = 0;
			for (int r = 0; r < ny; ++r){
				for (int c = 0; c < nx; ++c) {
					float last_val = rdata[nz/2*nxy+r*nx+c];
					for (int s = nz-1; s >= nz/2; --s) {
						idx = s*nxy+r*nx+c;
						tmp = rdata[idx];
						rdata[idx] = last_val;
						last_val = tmp;
					}
				}
			}
		}


	}
}

void AutoMaskAsymUnit::process_inplace(EMData* image) {
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	int ox = nx/2;
	int oy = ny/2;
	int oz = nz/2;

	Symmetry3D* sym = Factory<Symmetry3D>::get((string)params["sym"]);
	int au = params.set_default("au",0);

	float *d = image->get_data();
	for(int k = 0; k < nz; ++k ) {
		for(int j = 0; j < ny; ++j ) {
			for (int i = 0; i< nx; ++i, ++d) {
				//cout << i << " " << j << " " << k << endl;
				Vec3f v(i-ox,j-oy,k-oz);
// 				v.normalize();
				int a = sym->point_in_which_asym_unit(v);
				if (au == -1) {
					*d = (float)a;
				} else {
					if ( a == au ) *d = 1;
					else *d = 0;
				}
			}
		}
	}

	delete sym;

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
#ifdef	_WIN32
			if (_hypot(l - d_ny / 2, j - d_nx / 2) >= d_ny / 2) {
#else
			if (hypot(l - d_ny / 2, j - d_nx / 2) >= d_ny / 2) {
#endif
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

	Randnum * randnum = Randnum::Instance();
	if(params.has_key("seed")) {
		randnum->set_seed((int)params["seed"]);
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	image->ap2ri();
	float *rdata = image->get_data();

	size_t k = 0;
	float half_nz = 0;
	if (nz > 1) {
		half_nz = nz / 2.0f;
	}

	const float sqrt_2 = sqrt((float) 2);

	float r;
	for (int h = 0; h < nz; h++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i += 2, k += 2) {
				r = (Util::hypot3(i / 2.0f, j - ny / 2.0f, h - half_nz));
//				r = sqrt(Util::hypot3(i / 2.0f, j - ny / 2.0f, h - half_nz)); // I don't think this sqrt was supposed to be here --steve
				r = (r - x0) / dx;
				int l = 0;
				if (interpolation) {
					l = Util::fast_floor(r);
				}
				else {
					l = Util::fast_floor(r + 0.5f);
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
				f = randnum->get_gauss_rand(sqrt(f), sqrt(f) / 3);
				float a = randnum->get_frand(0.0f, (float)(2 * M_PI));
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
				size_t cur_z = z * nx * ny;

				for (int y = 1; y < ny - 1; y++) {
					size_t cur_y = y * nx + cur_z;

					for (int x = 1; x < nx - 1; x++) {
						size_t j = x + cur_y;

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

	int int_shift_only = params.set_default("int_shift_only",1);

	if ((float)image->get_attr("sigma")==0.0f) return;		// Can't center a constant valued image

	FloatPoint com = image->calc_center_of_mass();

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	if (int_shift_only) {
		int dx = -(int)(floor(com[0] + 0.5f) - nx / 2);
		int dy = -(int)(floor(com[1] + 0.5f) - ny / 2);
		int dz = 0;
		if (nz > 1) {
			dz = -(int)(floor(com[2] + 0.5f) - nz / 2);
		}
		image->translate(dx, dy, dz);

		Transform t;
		t.set_trans((float)dx,(float)dy,(float)dz);

		if (nz > 1) {
			image->set_attr("xform.align3d",&t);
		} else {
			image->set_attr("xform.align2d",&t);
		}
	}
	else {
		float dx = -(com[0] - nx / 2);
		float dy = -(com[1] - ny / 2);
		float dz = 0;
		if (nz > 1) {
			dz = -(com[2] - nz / 2);
		}
		image->translate(dx, dy, dz);

		Transform t;
		t.set_trans(dx,dy,dz);

		if (nz > 1) {
			image->set_attr("xform.align3d",&t);
		} else {
			image->set_attr("xform.align2d",&t);
		}
	}
}

void PhaseToMassCenterProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int int_shift_only = params.set_default("int_shift_only",1);

	vector<float> pcog = image->phase_cog();

	int dims = image->get_ndim();

	if (int_shift_only) {
		int dx=-int(pcog[0]+0.5f),dy=0,dz=0;
		if ( dims >= 2 ) dy = -int(pcog[1]+0.5);
		if ( dims == 3 ) dz = -int(pcog[2]+0.5);

		Transform t;
		t.set_trans((float)dx,(float)dy,(float)dz);
		if (dims == 3) image->set_attr("xform.align3d",&t);
		else if (dims == 2) image->set_attr("xform.align2d",&t);

		image->translate(dx,dy,dz);
	} else  {
		float dx=-pcog[0],dy=0.0,dz=0.0;
		if ( dims >= 2 ) dy = -pcog[1];
		if ( dims == 3 ) dz = -pcog[2];
		image->translate(dx,dy,dz);

		Transform t;
		t.set_trans(dx,dy,dz);
		if (dims == 3) image->set_attr("xform.align3d",&t);
		else if (dims == 2) image->set_attr("xform.align2d",&t);
	}
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
	EMData* aligned = image->align("translational", 0, params1);
	if ( image->get_ndim() == 3 ) {
		Transform* t = aligned->get_attr("xform.align3d");
		image->translate(t->get_trans());
		image->set_attr("xform.align3d",t);
		delete t;
	}
	else {
		// assumption is the image is 2D which may be  false
		Transform* t = aligned->get_attr("xform.align2d");
		image->translate(t->get_trans());
		image->set_attr("xform.align2d",t);
		delete t;
	}

	delete aligned;

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
		ctf = image_ctf->compute_1d(image->get_ysize(),1.0f/(image_ctf->apix*image->get_ysize()), Ctf::CTF_WIENER_FILTER, &sf);
	}
	else {
		ctf = image_ctf->compute_1d(image->get_ysize(),1.0f/(image_ctf->apix*image->get_ysize()), Ctf::CTF_SNR, &sf);
	}

	if(image_ctf) {delete image_ctf; image_ctf=0;}

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
	vector< vector< Transform > > transforms(sym_num);
	vector< float* > symvals(sym_num);
	for (int i =0; i < sym_num; i++) {
		vector<Transform> sym_transform =  Symmetry3D::get_symmetries(sym_list[i]);
		transforms[i] = sym_transform;
		symvals[i] = new float[sym_transform.size()]; // new float(nsym);
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
				size_t index = k * nx * ny + j * nx + i;
				float val = sdata[ index ];
				float bestmean = val, bestsymlevel = FLT_MAX;
				int bestsym = 0;
				for( int sym = 0; sym< sym_num; sym++) {
					int cur_sym_num = transforms[sym].size();
					float *symval = symvals[sym];
					// first find out all the symmetry related location values
					for( int s = 0; s < cur_sym_num; s++){
						Transform r = transforms[sym][s];
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

	float xo = image->get_attr("origin_x");
	float yo = image->get_attr("origin_y");
	float zo = image->get_attr("origin_z");

	float xom = msk->get_attr("origin_x");
	float yom = msk->get_attr("origin_y");
	float zom = msk->get_attr("origin_z");

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
					size_t idx = 0;
					float xc;
					int im;
					for (int i = 0; i < nx; i++) {
						xc = xo + i * apix;
						idx = k * nxy + j * nx + i;
						if (xc <= xom || xc >= xom + xm * apixm) {
							dp[idx] = 0;
						}
						else {
							im = (int) ((xc - xom) / apixm);
							if (dpm[km * xm * ym + jm * xm + im] <= 0) {
								dp[idx] = 0;
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

void MatchSFProcessor::create_radial_func(vector < float >&radial_mask,EMData *image) const {
	// The radial mask comes in with the existing radial image profile
	// The radial mask runs from 0 to the 1-D Nyquist (it leaves out the corners in Fourier space)

	EMData *to = params["to"];
	XYData *sf = new XYData();
	float apixto = to->get_attr("apix_x");


	if (to->is_complex()) {
		vector<float> rd=to->calc_radial_dist(to->get_ysize(),0,0.5,1);
		for (size_t i=0; i<rd.size(); i++) {
			sf->set_x(i,i/(apixto*2.0f*rd.size()));
			sf->set_y(i,rd[i]);
		}
	}
	else {
		EMData *tmp=to->do_fft();
		vector<float> rd=tmp->calc_radial_dist(to->get_ysize(),0,0.5,1);
		for (size_t i=0; i<rd.size(); i++) {
			sf->set_x(i,i/(apixto*2.0f*rd.size()));
			sf->set_y(i,rd[i]);
		}
		delete tmp;
	}

	float apix=image->get_attr("apix_x");

	int n = radial_mask.size();
	for (int i=0; i<n; i++) {
		if (radial_mask[i]>0) radial_mask[i]= sqrt(sf->get_yatx(i/(apix*2.0f*n))/radial_mask[i]);
	}

	delete sf;
}

void SetSFProcessor::create_radial_func(vector < float >&radial_mask,EMData *image) const {
	// The radial mask comes in with the existing radial image profile
	// The radial mask runs from 0 to the 1-D Nyquist (it leaves out the corners in Fourier space)

	XYData *sf = params["strucfac"];
	if(params.has_key("apix")) {
		image->set_attr("apix_x", (float)params["apix"]);
		image->set_attr("apix_y", (float)params["apix"]);
		image->set_attr("apix_z", (float)params["apix"]);
	}

	float apix=image->get_attr("apix_x");

	int n = radial_mask.size();
	for (int i=0; i<n; i++) {
		if (radial_mask[i]>0) radial_mask[i]= n*n*n*sqrt(sf->get_yatx(i/(apix*2.0f*n))/radial_mask[i]);
	}

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
	float r = 0.0f;
	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i, ++dat) {
				r =
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
	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i, ++dat) {
				r =
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
			size_t k2 = k * nxy;
			for (int j = 1; j < ny - 1; j++) {
				size_t l = j * nx + k2 + 1;

				for (int i = 1; i < nx - 1; i++) {
					if (dat[l] >= threshold || dat2[l]) {
						if (dat2[l - 1] || dat2[l + 1] ||
							dat2[l - nx] || dat2[l + nx] || dat2[l - nxy] || dat2[l + nxy]) {
							dat2[l] = 1.0f;
							done = false;
						}
					}
					++l;
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
	size_t idx;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			int j2 = j * nx + i;
			int k0 = 0;
			for (int k = 0; k < nz; k++) {
				idx = j2 + k * nxy;
				if (dat2[idx]) {
					k0 = k;
					break;
				}
			}

			if (k0 != nz) {
				int k1 = nz - 1;
				for (int k = nz - 1; k >= 0; k--) {
					idx = j2 + k * nxy;
					if (dat2[idx]) {
						k1 = k;
						break;
					}
				}

				for (int k = k0 + 1; k < k1; k++) {
					idx = j2 + k * nxy;
					dat2[idx] = 1.0f;
				}
			}
		}
	}

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < nz; j++) {
			size_t j2 = j * nxy + i;
			int k0 = 0;
			for (int k = 0; k < ny; k++) {
				idx = k * nx + j2;
				if (dat2[idx]) {
					k0 = k;
					break;
				}
			}

			if (k0 != ny) {
				int k1 = ny - 1;
				for (int k = ny - 1; k >= 0; k--) {
					idx = k * nx + j2;
					if (dat2[idx]) {
						k1 = k;
						break;
					}
				}

				for (int k = k0 + 1; k < k1; k++) {
					idx = k * nx + j2;
					dat2[idx] = 1.0f;
				}
			}
		}
	}

	for (int i = 0; i < ny; i++) {
		for (int j = 0; j < nz; j++) {
			size_t j2 = i * nx + j * nxy;
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

	size_t l = 0;
	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i) {
				if (dat[l] > t) {
					dat2[l] = 1.0f;
				}
				++l;
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

	for (int k = 1; k < nz - 1; ++k) {
		for (int j = 1; j < ny - 1; ++j) {
			size_t l = j * nx + k * nxy + 1;
			for (int i = 1; i < nx - 1; ++i, ++l) {
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

	if (image->get_ndim() != 3) {
		throw ImageDimensionException("This processor was only ever designed to work on 3D images.");
	}

	/*
	 The mask writing functionality was removed to comply with an EMAN2 policy which dictates that file io is not allowed from within a processor
	 To get around this just use the return_mask parameter.
	string mask_output = params.set_default("write_mask", "");
	if ( mask_output != "") {
		if (Util::is_file_exist(mask_output) ) throw InvalidParameterException("The mask output file name already exists. Please remove it if you don't need it.");
		if (!EMUtil::is_valid_filename(mask_output)) throw InvalidParameterException("The mask output file name type is invalid or unrecognized");
	}
	*/

	int radius=0;
	if (params.has_key("radius")) {
		radius = params["radius"];
	}
	int nmaxseed=0;
	if (params.has_key("nmaxseed")) {
		nmaxseed = params["nmaxseed"];
	}
	float threshold = params["threshold"];
	int nshells = params["nshells"];
	int nshellsgauss = params["nshellsgauss"];
	int verbose=params.set_default("verbose",0);

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	int nxy=nx*ny;

	EMData *amask = new EMData();
	amask->set_size(nx, ny, nz);

	float *dat = image->get_data();
	float *dat2 = amask->get_data();
	int i,j,k;
	size_t l = 0;

	// Seeds with the highest valued pixels
	if (nmaxseed>0) {
		vector<Pixel> maxs=image->calc_n_highest_locations(nmaxseed);
		
		for (vector<Pixel>::iterator i=maxs.begin(); i<maxs.end(); i++) {
			amask->set_value_at((*i).x,(*i).y,(*i).z,1.0);
			if (verbose) printf("Seed at %d,%d,%d (%1.3f)\n",(*i).x,(*i).y,(*i).z,(*i).value);
		}
	}
	
	// Seeds with a sphere
	if (radius>0) {
		// start with an initial sphere
		for (k = -nz / 2; k < nz / 2; ++k) {
			for (j = -ny / 2; j < ny / 2; ++j) {
				for (i = -nx / 2; i < nx / 2; ++i,++l) {
					if (abs(k) > radius || abs(j) > radius || abs(i) > radius) continue;
					if ( (k * k + j * j + i * i) > (radius*radius) || dat[l] < threshold) continue;
					dat2[l] = 1.0f;
				}
			}
		}
	}

	
	// iteratively 'flood fills' the map... recursion would be better
	int done=0;
	int iter=0;
	while (!done) {
		iter++;
		done=1;
		if (verbose && iter%10==0) printf("%d iterations\n",iter);
		for (k=1; k<nz-1; ++k) {
			for (j=1; j<ny-1; ++j) {
				for (i=1; i<nx-1; ++i) {
					l=i+j*nx+k*nx*ny;
					if (dat2[l]) continue;
					if (dat[l]>threshold && (dat2[l-1]||dat2[l+1]||dat2[l+nx]||dat2[l-nx]||dat2[l-nxy]||dat2[l+nxy])) {
						dat2[l]=1.0;
						done=0;
					}
				}
			}
		}
	}

	amask->update();

	if (verbose) printf("extending mask\n");
	amask->process_inplace("mask.addshells.gauss", Dict("val1", nshells, "val2", nshellsgauss));

	bool return_mask = params.set_default("return_mask",false);
	if (return_mask) {
		// Yes there is probably a much more efficient way of getting the mask itself, but I am only providing a stop gap at the moment.
		memcpy(dat,dat2,image->get_size()*sizeof(float));
	} else {
		image->mult(*amask);
	}

	// EMAN2 policy is not to allow file io from with a processor
	//if (mask_output != "") {
	//	amask->write_image(mask_output);
	//}


	delete amask;
}

void IterBinMaskProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	float val1 = params["val1"];
	float val2 = params["val2"];

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	EMData *image2 = new EMData(nx,ny,nz);

	// Got a compile warning complaining that these things were never used. Hope I didn't break anything - apologies if so (d.woolford)
//	float *dat = image->get_data();
//	float *dat2 = image2->get_data();


	float *d = image->get_data();
	float *d2 = image2->get_data();

	const int nxy = nx * ny;
	size_t size = nx * ny * nz;

	// TODO: THIS IS EXTREMELY INEFFICIENT
	if (nz != 1) {
		for (int l = 1; l <= (int) val1+val2; ++l) {
			for (size_t i=0; i<size; i++) d2[i]=d[i];
			for (int k = 1; k < nz - 1; ++k) {
				for (int j = 1; j < ny - 1; ++j) {
					for (int i = 1; i < nx - 1; ++i) {
						size_t t = i + j*nx+k*nx*ny;
						if (d[t]) continue;
						if (d2[t - 1] || d2[t + 1] || d2[t + nx] || d2[t - nx] || d2[t + nxy] || d2[t - nxy]) d[t] = (float) l + 1;
					}
				}
			}
		}
	}
	else {
		for (int l = 1; l <= (int) val1+val2; ++l) {
			for (size_t i=0; i<size; i++) d2[i]=d[i];
			for (int j = 1; j < ny - 1; ++j) {
				for (int i = 1; i < nx - 1; ++i) {
					size_t t = i + j * nx;
					if (d[t]) continue;
					if (d2[t - 1] || d2[t + 1] || d2[t + nx] || d2[t - nx]) d[t] = (float) l + 1;
				}
			}
		}
	}

	vector<float> vec;
	for (int i=0; i<val1+2; i++) vec.push_back(1.0);
	for (int i=0; i<val2; i++) {
		vec.push_back(exp(-pow(2.0f*i/(val2),2.0f)));
//		printf("%f\n",exp(-pow(2.0*i/(val1-val2),2.0)));
	}
	for (size_t i = 0; i < size; i++) if (d[i]) d[i]=vec[(int)d[i]];

	image->update();
	delete image2;
}

EMData* DirectionalSumProcessor::process(const EMData* const image ) {
	string dir = params.set_default("direction", "");
	if ( dir == "" || ( dir != "x" && dir != "y" && dir != "z" ) )
		throw InvalidParameterException("The direction parameter must be either x, y, or z");

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	// compress one of the dimensions
	if ( dir == "x" ) nx = 1;
	else if ( dir == "y" ) ny = 1;
	else if ( dir == "z" ) nz = 1;

	EMData* ret = new EMData;
	ret->set_size(nx,ny,nz);
	ret->to_zero();

	float* d = image->get_data();
	for(int k = 0; k < image->get_zsize(); ++k ) {
		for(int j = 0; j < image->get_ysize(); ++j ) {
			for(int i = 0; i < image->get_xsize(); ++i, ++d ) {
				if ( dir == "x" ) {
					float v = ret->get_value_at(0,j,k);
					ret->set_value_at(0,j,k,*d+v);
				}else if ( dir == "y" ) {
					float v = ret->get_value_at(i,0,k);
					ret->set_value_at(i,0,k,*d+v);
				}
				else if ( dir == "z" ) {
					float v = ret->get_value_at(i,j,0);
					ret->set_value_at(i,j,0,*d+v);
				}
			}
		}
	}
	ret->update();
	return ret;
}

void TestImageProcessor::preprocess(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	nx = image->get_xsize();
	ny = image->get_ysize();
	nz = image->get_zsize();
}


void TestImageFourierNoiseGaussian::process_inplace(EMData* image)
{
	if (!image->is_complex()) {
		int nx = image->get_xsize();
		int offset = 2 - nx%2;

		image->set_size(nx+offset,image->get_ysize(),image->get_zsize());
		image->set_complex(true);
		if (1 == offset) image->set_fftodd(true);
		else image->set_fftodd(false);
		image->set_fftpad(true);
	}
	image->ri2ap();

	float sigma = params.set_default("sigma",.25f);

	float * d = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nxy = image->get_ysize()*nx;
	int nzon2 = image->get_zsize()/2;
	int nyon2 = image->get_ysize()/2;
	float rx, ry, rz, length, amp, phase;
	int twox;
	for (int z = 0; z< image->get_zsize(); ++z) {
		for (int y = 0; y < image->get_ysize(); ++y) {
			for (int x = 0; x < image->get_xsize()/2; ++x) {
				rx = (float)x;
				ry = (float)nyon2 - (float)y;
				rz = (float)nzon2 - (float)z;
				length = sqrt(rx*rx + ry*ry + rz*rz);
				amp = exp(-sigma*length);
				phase = Util::get_frand(0,1)*2*M_PI;

				twox = 2*x;
				size_t idx1 = twox + y*nx+z*nxy;
				size_t idx2 = idx1 + 1;
				d[idx1] = amp;
				d[idx2] = phase;

			}
		}
	}

	image->ap2ri();
	if (image->get_ndim() == 2) {
		bool yodd = image->get_ysize() % 2 == 1;

		int yit = image->get_ysize()/2-1;
		int offset = 1;
		if (yodd) {
			offset = 0;
		}
		for (int y = 0; y < yit; ++y) {
			int bot_idx = (y+offset)*nx;
			int top_idx = (ny-1-y)*nx;
			float r1 = d[bot_idx];
			float i1 = d[bot_idx+1];
			float r2 = d[top_idx];
			float i2 = d[top_idx+1];
			float r = (r1 + r2)/2.0f;
			float i = (i1 + i2)/2.0f;
			d[bot_idx] = r;
			d[top_idx] = r;
			d[bot_idx+1] = i;
			d[top_idx+1] = -i;

			bot_idx = (y+offset)*nx+nx-2;
			top_idx = (ny-1-y)*nx+nx-2;
			r1 = d[bot_idx];
			i1 = d[bot_idx+1];
			r2 = d[top_idx];
			i2 = d[top_idx+1];
			r = (r1 + r2)/2.0f;
			i = (i1 + i2)/2.0f;
			d[bot_idx] = r;
			d[top_idx] = r;
			d[bot_idx+1] = i;
			d[top_idx+1] = -i;
		}

		d[1] = 0; // 0 phase for this componenet
		d[nx-1] = 0; // 0 phase for this component
		d[ny/2*nx+nx-1] = 0;// 0 phase for this component
		d[ny/2*nx+1] = 0;// 0 phase for this component
	}

	if (image->get_ndim() != 1) image->process_inplace("xform.fourierorigin.tocorner");
	image->do_ift_inplace();
	image->depad();
}

#include <iostream>
using std::ostream_iterator;

void CTFSNRWeightProcessor::process_inplace(EMData* image) {
	if (params.has_key("noise")==false) throw InvalidParameterException("You must supply the noise argument");
	if (params.has_key("snr")==false) throw InvalidParameterException("You must supply the snr argument");

	float boost = params.set_default("boost",1.0f);

	if (!image->is_complex()) {
		image->do_fft_inplace();
	}
	EMData* cpy = image->copy();
	cpy->ri2inten();
	vector<float> sf = cpy->calc_radial_dist(cpy->get_ysize()/2,0.0,1.0,1);
	transform(sf.begin(),sf.end(),sf.begin(),sqrtf);
	delete cpy;

	image->ri2ap();

	vector<float> noise = params["noise"];
	vector<float> snr = params["snr"];

// 	copy(snr.begin(), snr.end(), ostream_iterator<float>(cout, "\n"));
// 	copy(noise.begin(), noise.end(), ostream_iterator<float>(cout, "\n"));

	for(vector<float>::iterator it = noise.begin(); it != noise.end(); ++it){
		if ((*it) == 0) *it = 1;
	}
	for(vector<float>::iterator it = snr.begin(); it != snr.end(); ++it){
		if ((*it) < 0) *it = 0;
	}
	// Subtract the mean from the data and store it in data_mm
	transform(snr.begin(),snr.end(),noise.begin(),snr.begin(),std::multiplies<float>());
	transform(snr.begin(),snr.end(),snr.begin(),sqrtf);
// 	copy(snr.begin(), snr.end(), ostream_iterator<float>(cout, "\n"));
// 	copy(noise.begin(), noise.end(), ostream_iterator<float>(cout, "\n"));

	int i = static_cast<int>(snr.size());

	float * d = image->get_data();
	int nx = image->get_xsize();
// 	int ny = image->get_ysize();
	int nxy = image->get_ysize()*nx;
	int nzon2 = image->get_zsize()/2;
	int nyon2 = image->get_ysize()/2;
	float rx, ry, rz, amp;
	int length;
	int twox;
	image->process_inplace("xform.fourierorigin.tocenter");
	for (int z = 0; z< image->get_zsize(); ++z) {
		for (int y = 0; y < image->get_ysize(); ++y) {
			for (int x = 0; x < image->get_xsize()/2; ++x) {
				rx = (float)x;
				ry = (float)nyon2 - (float)y;
				rz = (float)nzon2 - (float)z;
				length = static_cast<int>(sqrt(rx*rx + ry*ry + rz*rz));

				twox = 2*x;
				size_t idx1 = twox + y*nx+z*nxy;
				if (length >= i || length >= (int)sf.size()) {
					d[idx1] = 0;
					continue;
				} else {
					amp = boost*snr[length];
// 					if (amp > 0) amp =sqrtf(amp);
// 					else amp = 0;
				}

				if (sf[length] == 0) {
					d[idx1] = 0;
					continue;
				}

// 				size_t idx2 = idx1 + 1;
// 				cout << d[idx1] << " " << sf[length] << endl;
				d[idx1] /= sf[length];
				if (d[idx1] < 0) {
					d[idx1] *= amp;
				}else {
					d[idx1] *= -amp;
				}
// 				d[idx2] = phase;

			}
		}
	}

	image->ap2ri();
	if (image->get_ndim() != 1) image->process_inplace("xform.fourierorigin.tocorner");
	image->do_ift_inplace();
	image->depad();
}

void TestImageFourierNoiseProfile::process_inplace(EMData * image) {

	if (params.has_key("profile")==false) throw InvalidParameterException("You must supply the profile argument");

	if (!image->is_complex()) {
		int nx = image->get_xsize();
		int offset = 2 - nx%2;

		image->set_size(nx+offset,image->get_ysize(),image->get_zsize());
		image->set_complex(true);
		if (1 == offset) image->set_fftodd(true);
		else image->set_fftodd(false);
		image->set_fftpad(true);
	}
	image->to_zero();
	image->ri2ap();

	vector<float> profile = params["profile"];
	transform(profile.begin(),profile.end(),profile.begin(),sqrtf);

	int i = static_cast<int>(profile.size());

	float * d = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nxy = image->get_ysize()*nx;
	int nzon2 = image->get_zsize()/2;
	int nyon2 = image->get_ysize()/2;
	float rx, ry, rz, amp, phase;
	int length;
	int twox;
	for (int z = 0; z< image->get_zsize(); ++z) {
		for (int y = 0; y < image->get_ysize(); ++y) {
			for (int x = 0; x < image->get_xsize()/2; ++x) {
				rx = (float)x;
				ry = (float)nyon2 - (float)y;
				rz = (float)nzon2 - (float)z;
				length = static_cast<int>(sqrt(rx*rx + ry*ry + rz*rz));

				twox = 2*x;
				size_t idx1 = twox + y*nx+z*nxy;
				size_t idx2 = idx1 + 1;


				if (length >= i) {
					d[idx1] = 0;
					d[idx2] = 0;
					continue;
				}
				amp = profile[length];
				phase = Util::get_frand(0,1)*2*M_PI;


				d[idx1] = amp;
				d[idx2] = phase;

			}
		}
	}

	image->ap2ri();
	if (image->get_ndim() == 2) {
		bool yodd = image->get_ysize() % 2 == 1;

		int yit = image->get_ysize()/2-1;
		int offset = 1;
		if (yodd) {
			offset = 0;
		}
		for (int y = 0; y < yit; ++y) {
			int bot_idx = (y+offset)*nx;
			int top_idx = (ny-1-y)*nx;
			float r1 = d[bot_idx];
			float i1 = d[bot_idx+1];
			float r2 = d[top_idx];
			float i2 = d[top_idx+1];
			float r = (r1 + r2)/2.0f;
			float i = (i1 + i2)/2.0f;
			d[bot_idx] = r;
			d[top_idx] = r;
			d[bot_idx+1] = i;
			d[top_idx+1] = -i;

			bot_idx = (y+offset)*nx+nx-2;
			top_idx = (ny-1-y)*nx+nx-2;
			r1 = d[bot_idx];
			i1 = d[bot_idx+1];
			r2 = d[top_idx];
			i2 = d[top_idx+1];
			r = (r1 + r2)/2.0f;
			i = (i1 + i2)/2.0f;
			d[bot_idx] = r;
			d[top_idx] = r;
			d[bot_idx+1] = i;
			d[top_idx+1] = -i;
		}

		d[1] = 0; // 0 phase for this componenet
		d[nx-1] = 0; // 0 phase for this component
		d[ny/2*nx+nx-1] = 0;// 0 phase for this component
		d[ny/2*nx+1] = 0;// 0 phase for this component
	}

	if (image->get_ndim() != 1) image->process_inplace("xform.fourierorigin.tocorner");
	image->do_ift_inplace();
	image->depad();
}

void TestImageLineWave::process_inplace(EMData * image)
{
	preprocess(image);

	float period = params.set_default("period",10.0f);
	int n = image->get_xsize()*image->get_ysize()*image->get_zsize();

	for(int i = 0; i < n; ++i) {
		float x = fmod((float)i,period);
		x /= period;
		x = (float)sin(x*EMConsts::pi*2.0);
		image->set_value_at_fast(i,x);
	}
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
	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i, ++dat) {
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

	float m = params.set_default("m", 1.0f);
	float b = params.set_default("b", 0.0f);

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
		throw ImageDimensionException("works for 2D images only");
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	image->to_zero();

	for (int i=0; i<100; i++) {
		int x=static_cast<int>( nx/2+nx/6.0*sin(i*2.0*3.14159/100.0) );
		int y=ny/4+i*ny/200;
		for (int xx=x-nx/10; xx<x+nx/10; xx++) {
			for (int yy=y-ny/10; yy<y+ny/10; yy++) {
#ifdef	_WIN32
				(*image)(xx,yy)+=exp(-pow(static_cast<float>(_hypot(xx-x,yy-y))*30.0f/nx,2.0f))*(sin(static_cast<float>((xx-x)*(yy-y)))+.5f);
#else
				(*image)(xx,yy)+=exp(-pow(static_cast<float>(hypot(xx-x,yy-y))*30.0f/nx,2.0f))*(sin(static_cast<float>((xx-x)*(yy-y)))+.5f);
#endif
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
	float norm  = 1.0f/ ( x_sigma*sr2pi );
	if (ny > 1) {
		norm *= 1.0f/ ( y_sigma*sr2pi );
		if (nz > 1) norm  *= 1.0f/ ( z_sigma*sr2pi );
	}

	float z, y, x, sum, val;
	for (int iz=0; iz < nz; ++iz) {
		z = static_cast<float>(iz) - z_center;
		for (int iy=0; iy < ny; ++iy) {
			y = static_cast<float>(iy) - y_center;
			for (int ix=0; ix < nx; ++ix) {
				x = static_cast<float>(ix) - x_center;
				sum = x*x/x_twosig2 + y*y/y_twosig2 + z*z/z_twosig2;
				val = norm*exp(-sum);
				(*image)(ix,iy,iz) = val;
			}
		}
	}
	image->update();
}

void TestImageSphericalWave::process_inplace(EMData * image)
{
	preprocess(image);

	if(!params.has_key("wavelength")) {
		LOGERR("%s wavelength is required parameter", get_name().c_str());
		throw InvalidParameterException("wavelength parameter is required.");
	}
	float wavelength = params["wavelength"];

	float phase = 0;
	if(params.has_key("phase")) {
		phase = params["phase"];
	}

	float x = (float)(nx/2);
	if (params.has_key("x")) x=params["x"];
	float y = (float)(ny/2);
	if (params.has_key("y")) y=params["y"];
	float z = (float)(nz/2);
	if (params.has_key("z")) z=params["z"];

	int ndim = image->get_ndim();

	if(ndim==2) {	//2D
		for(int j=0; j<ny; ++j) {
			for(int i=0; i<nx; ++i) {
#ifdef _WIN32
				float r=_hypotf(x-(float)i,y-(float)j);
#else
				float r=hypot(x-(float)i,y-(float)j);
#endif	//_WIN32
				if (r<.5) continue;
				image->set_value_at(i,j,cos(2*(float)pi*r/wavelength+phase)/r);
			}
		}
	}
	else {	//3D
		for(int k=0; k<nz; ++k) {
			for(int j=0; j<ny; ++j) {
				for(int i=0; i<nx; ++i) {
					float r=Util::hypot3(x-(float)i,y-(float)j,z-(float)k);
					if (r<.5) continue;
					image->set_value_at(i,j,k,cos(2*(float)pi*r/wavelength+phase)/(r*r));
				}
			}
		}
	}

	image->update();
}


void TestImageSinewave::process_inplace(EMData * image)
{
	preprocess(image);

	if(!params.has_key("wavelength")) {
		LOGERR("%s wavelength is required parameter", get_name().c_str());
		throw InvalidParameterException("wavelength parameter is required.");
	}
	float wavelength = params["wavelength"];

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
			*dat = sin(i*(2.0f*M_PI/wavelength) + phase);
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
					*dat = sin((i*sin((180-alpha)*M_PI/180)+j*cos((180-alpha)*M_PI/180))*(2.0f*M_PI/wavelength) + phase);
				}
				else if(axis.compare("y")==0 || axis.compare("Y")==0) {
					*dat = sin(j*(2.0f*M_PI/wavelength) + phase);
				}
				else {
					*dat = sin(i*(2.0f*M_PI/wavelength) + phase);
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
						*dat = sin(k*(2.0f*M_PI/wavelength) + phase);
					}
					else if(axis.compare("y")==0 || axis.compare("Y")==0) {
						*dat = sin(j*(2.0f*M_PI/wavelength) + phase);
					}
					else {
						*dat = sin(i*(2.0f*M_PI/wavelength) + phase);
					}
				}
			}
		}

		if(az != 0 || alt != 0 || phi != 0) {
			Dict d("type","eman");
			d["az"] = az; d["phi"] = phi; d["alt"] = alt;
			image->transform(Transform(d));
		}
	}

	image->update();
}

void TestImageSinewaveCircular::process_inplace(EMData * image)
{
	preprocess(image);

	float wavelength = params["wavelength"];
	string axis = (const char*)params["axis"];
	float c = params["c"];
	float phase = params["phase"];

	float *dat = image->get_data();
	float r; //this is the distance of pixel from the image center(nx/2, ny/2, nz/2)
	float x2, y2, z2; //this is the coordinates of this pixel from image center
	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i, ++dat) {
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
				*dat = sin( r * (2.0f*M_PI/wavelength) - phase*180/M_PI);
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
	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i, ++dat) {
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

	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i, ++dat) {
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

void TestImageHollowEllipse::process_inplace(EMData * image)
{
	preprocess(image);

	float width = params.set_default("width",2.0f);

	float a2 = params.set_default("a",nx/2.0f-1.0f);
	float b2 = params.set_default("b",ny/2.0f-1.0f);
	float c2 = params.set_default("c",nz/2.0f-1.0f);

	float a1 = params.set_default("xwidth",a2-width);
	float b1 = params.set_default("ywidth",b2-width);
	float c1 = params.set_default("zwidth",c2-width);

	float fill = params.set_default("fill",1.0f);
	Transform* t;
	if (params.has_key("transform")) {
		t = params["transform"];
	} else {
		t = new Transform;
	}


	int mz = 2*(int)c2+1;
	if ( nz < mz ) mz = nz;
	int my = 2*(int)b2+1;
	if ( ny < my ) my = ny;
	int mx = 2*(int)a2+1;
	if ( nx < mx ) mx = nx;

	float ai1 = 1/(a1*a1);
	float bi1 = 1/(b1*b1);
	float ci1 = 1/(c1*c1);

	float ai2 = 1/(a2*a2);
	float bi2 = 1/(b2*b2);
	float ci2 = 1/(c2*c2);

	Vec3f origin(nx/2,ny/2,nz/2);

	float x2, y2, z2, r1,r2;
	int xl, yl, zl;
	for (int k = 0; k < mz; ++k) {
		for (int j = 0; j < my; ++j) {
			for (int i = 0; i < mx; ++i) {
				x2 = (float)(i - mx/2);
				y2 = (float)(j - my/2);
				z2 = (float)(k - mz/2);
				r1 = (x2*x2)*ai1 + (y2*y2)*bi1 + (z2*z2)*ci1;
				r2 = (x2*x2)*ai2 + (y2*y2)*bi2 + (z2*z2)*ci2;
				if (r2 <= 1 && r1 >= 1) {

					if (t != 0) {
						Vec3f v(x2,y2,z2);
						v = (*t)*v;
						v += origin;

						// THIS ISN'T THE BEST STRATEGY BUT IT'S A STOP GAP. A FLOOD FILL IS PROBABLY BETTER
						// I fill in 3x3 cubes to make sure there are no gaps...

						for( int kk = -1; kk <= 1; ++kk)
							for( int jj = -1; jj <= 1; ++jj)
								for( int ii = -1; ii <= 1; ++ii) {
									xl = (int)v[0]+ii;
									yl = (int)v[1]+jj;
									zl = (int)v[2]+kk;
									if (xl >= 0 && xl < nx && yl >= 0 && yl < ny && zl >= 0 && zl < nz)
										image->set_value_at(xl,yl,zl,1.0);
								}
					} else {
						image->set_value_at((int)x2+nx/2,(int)y2+ny/2,(int)z2+nz/2,fill);
					}
				}
			}
		}
	}

	delete t;

	image->update();
}

void TestImageEllipse::process_inplace(EMData * image)
{
	preprocess(image);


	float a = params.set_default("a",nx/2.0f-1.0f);
	float b = params.set_default("b",ny/2.0f-1.0f);
	float c = params.set_default("c",nz/2.0f-1.0f);
	float fill = params.set_default("fill",1.0f);
	//bool hollow = params.set_default("hollow",false);
	Transform* t;
	if (params.has_key("transform")) {
		t = params["transform"];
	} else {
		t = new Transform;
	}


	int mz = 2*(int)c+1;
	if ( nz < mz ) mz = nz;
	int my = 2*(int)b+1;
	if ( ny < my ) my = ny;
	int mx = 2*(int)a+1;
	if ( nx < mx ) mx = nx;


	float ai = 1/(a*a);
	float bi = 1/(b*b);
	float ci = 1/(c*c);

	Vec3f origin(nx/2,ny/2,nz/2);

	float x2, y2, z2, r;
	int xl, yl, zl;
	for (int k = 0; k < mz; ++k) {
		for (int j = 0; j < my; ++j) {
			for (int i = 0; i < mx; ++i) {
				x2 = (float)(i - mx/2);
				y2 = (float)(j - my/2);
				z2 = (float)(k - mz/2);
				r = (x2*x2)*ai + (y2*y2)*bi + (z2*z2)*ci;
				if (r <= 1) {

					if (t != 0) {
						Vec3f v(x2,y2,z2);
						v = (*t)*v;
						v += origin;

						// THIS ISN'T THE BEST STRATEGY BUT IT'S A STOP GAP. A FLOOD FILL IS PROBABLY BETTER
						// I fill in 3x3 cubes to make sure there are no gaps...

						for( int kk = -1; kk <= 1; ++kk)
							for( int jj = -1; jj <= 1; ++jj)
								for( int ii = -1; ii <= 1; ++ii) {
									xl = (int)v[0]+ii;
									yl = (int)v[1]+jj;
									zl = (int)v[2]+kk;
									if (xl >= 0 && xl < nx && yl >= 0 && yl < ny && zl >= 0 && zl < nz)
										image->set_value_at(xl,yl,zl,fill);
								}
					} else {
						image->set_value_at((int)x2+nx/2,(int)y2+ny/2,(int)z2+nz/2,fill);
					}
				}
			}
		}
	}

	delete t;

	image->update();
}

void TestImageNoiseUniformRand::process_inplace(EMData * image)
{
	preprocess(image);

	Randnum * r = Randnum::Instance();
	if(params.has_key("seed")) {
		r->set_seed((int)params["seed"]);
	}

	float *dat = image->get_data();
	size_t size = nx*ny*nz;
	for (size_t i=0; i<size; ++i) {
		dat[i] = r->get_frand();
	}

	image->update();
}

void TestImageNoiseGauss::process_inplace(EMData * image)
{
	preprocess(image);

	float sigma = params["sigma"];
	if (sigma<=0) { sigma = 1.0; }
	float mean = params["mean"];

	Randnum * r = Randnum::Instance();
	if (params.has_key("seed")) {
		r->set_seed((int)params["seed"]);
	}

	float *dat = image->get_data();
	size_t size = nx*ny*nz;
	for (size_t i=0; i<size; ++i) {
		dat[i] = r->get_gauss_rand(mean, sigma);
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
	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i, ++dat) {
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
//	float l=log((float)nx)/log(2.0f);
//	if (l!=floor(l)) throw ImageDimensionException("Wavelet transform size must be power of 2");
	if( !Util::IsPower2(nx) )  throw ImageDimensionException("Wavelet transform size must be power of 2");

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

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	size_t nxy = nx*ny;

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
 		 size_t beg     = 0;
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
			size_t beg     = 0;
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


float* TransformProcessor::transform(const EMData* const image, const Transform& t) const {

	ENTERFUNC;

	Transform inv = t.inverse();
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	int nxy = nx*ny;

	const float * const src_data = image->get_const_data();
	float *des_data = (float *) EMUtil::em_malloc(nx*ny*nz* sizeof(float));

	if (nz == 1) {
		Vec2f offset(nx/2,ny/2);
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				Vec2f coord(i-nx/2,j-ny/2);
				Vec2f soln = inv*coord;
				soln += offset;

				float x2 = soln[0];
				float y2 = soln[1];

				if (x2 < 0 || x2 >= nx || y2 < 0 || y2 >= ny ) {
					des_data[i + j * nx] = 0; // It may be tempting to set this value to the
					// mean but in fact this is not a good thing to do. Talk to S.Ludtke about it.
				}
				else {
					int ii = Util::fast_floor(x2);
					int jj = Util::fast_floor(y2);
					int k0 = ii + jj * nx;
					int k1 = k0 + 1;
					int k2 = k0 + nx;
					int k3 = k0 + nx + 1;

					if (ii == nx - 1) {
						k1--;
						k3--;
					}
					if (jj == ny - 1) {
						k2 -= nx;
						k3 -= nx;
					}

					float t = x2 - ii;
					float u = y2 - jj;

					des_data[i + j * nx] = Util::bilinear_interpolate(src_data[k0],src_data[k1], src_data[k2], src_data[k3],t,u);
				}
			}
		}
	}
	else {
		size_t l=0, ii, k0, k1, k2, k3, k4, k5, k6, k7;
		Vec3f offset(nx/2,ny/2,nz/2);
		float x2, y2, z2, tuvx, tuvy, tuvz;
		int ix, iy, iz;
		for (int k = 0; k < nz; ++k) {
			for (int j = 0; j < ny; ++j) {
				for (int i = 0; i < nx; ++i,++l) {
					Vec3f coord(i-nx/2,j-ny/2,k-nz/2);
					Vec3f soln = inv*coord;
					soln += offset;

					x2 = soln[0];
					y2 = soln[1];
					z2 = soln[2];

					if (x2 < 0 || y2 < 0 || z2 < 0 || x2 >= nx  || y2 >= ny  || z2>= nz ) {
						des_data[l] = 0;
					}
					else {
						ix = Util::fast_floor(x2);
						iy = Util::fast_floor(y2);
						iz = Util::fast_floor(z2);
						tuvx = x2-ix;
						tuvy = y2-iy;
						tuvz = z2-iz;
						ii = ix + iy * nx + iz * nxy;

						k0 = ii;
						k1 = k0 + 1;
						k2 = k0 + nx;
						k3 = k0 + nx+1;
						k4 = k0 + nxy;
						k5 = k1 + nxy;
						k6 = k2 + nxy;
						k7 = k3 + nxy;

						if (ix == nx - 1) {
							k1--;
							k3--;
							k5--;
							k7--;
						}
						if (iy == ny - 1) {
							k2 -= nx;
							k3 -= nx;
							k6 -= nx;
							k7 -= nx;
						}
						if (iz == nz - 1) {
							k4 -= nxy;
							k5 -= nxy;
							k6 -= nxy;
							k7 -= nxy;
						}

						des_data[l] = Util::trilinear_interpolate(src_data[k0],
								src_data[k1], src_data[k2], src_data[k3], src_data[k4],
								src_data[k5], src_data[k6],	src_data[k7], tuvx, tuvy, tuvz);
					}
				}
			}
		}
	}

	EXITFUNC;
	return des_data;
}

void TransformProcessor::assert_valid_aspect(const EMData* const image) const {
	int ndim = image->get_ndim();
	if (ndim != 2 && ndim != 3) throw ImageDimensionException("Transforming an EMData only works if it's 2D or 3D");

	if (! params.has_key("transform") ) throw InvalidParameterException("You must specify a Transform in order to perform this operation");
}

//void TransformProcessor::update_emdata_attributes(EMData* const p, const Dict& attr_dict, const float& scale) const {
//
//	float inv_scale = 1.0f/scale;
//
//	p->scale_pixel(1.0f/scale);
//
//	// According to Baker the origin attributes remain unchanged
////	vector<string> inv_scale_attrs;
////	inv_scale_attrs.push_back("origin_x");
////	inv_scale_attrs.push_back("origin_y");
////	inv_scale_attrs.push_back("origin_z");
////
////	for(vector<string>::const_iterator it = inv_scale_attrs.begin(); it != inv_scale_attrs.end(); ++it) {
////		if (attr_dict.has_key(*it)) {
////			p->set_attr(*it,(float) attr_dict[*it] * scale);
////		}
////	}
//
//	vector<string> scale_attrs;
//	scale_attrs.push_back("apix_x");
//	scale_attrs.push_back("apix_y");
//	scale_attrs.push_back("apix_z");
//
//
//	for(vector<string>::const_iterator it = scale_attrs.begin(); it != scale_attrs.end(); ++it) {
//		if (attr_dict.has_key(*it)) {
//			p->set_attr(*it,(float) attr_dict[*it] * inv_scale);
//		}
//	}
//
//}

EMData* TransformProcessor::process(const EMData* const image) {
	ENTERFUNC;

	assert_valid_aspect(image);

	Transform* t = params["transform"];

	EMData* p  = 0;
#ifdef EMAN2_USING_CUDA
	if (image->gpu_operation_preferred()) {
// 		cout << "cuda transform" << endl;
		float * m = new float[12];
		Transform inv = t->inverse();
		inv.copy_matrix_into_array(m);
		image->bind_cuda_texture();
		EMDataForCuda* tmp = emdata_transform_cuda(m,image->get_xsize(),image->get_ysize(),image->get_zsize());
		image->unbind_cuda_texture();
		delete [] m;
		if (tmp == 0) throw;

		p = new EMData();
		p->set_gpu_rw_data(tmp->data,tmp->nx,tmp->ny,tmp->nz);
		free(tmp);
	}
#endif
	if ( p == 0 ) {
		float* des_data = transform(image,*t);
		p = new EMData(des_data,image->get_xsize(),image->get_ysize(),image->get_zsize(),image->get_attr_dict());
	}

	// 	all_translation += transform.get_trans();

	float scale = t->get_scale();
	if (scale != 1.0) {
		p->scale_pixel(1.0f/scale);
//		update_emdata_attributes(p,image->get_attr_dict(),scale);
	}

	if(t) {delete t; t=0;}
	EXITFUNC;
	return p;
}

void TransformProcessor::process_inplace(EMData* image) {
	ENTERFUNC;

	assert_valid_aspect(image);

	Transform* t = params["transform"];

	// 	all_translation += transform.get_trans();
	bool use_cpu = true;
#ifdef EMAN2_USING_CUDA
	if (image->gpu_operation_preferred()) {
		use_cpu = false;
		float * m = new float[12];
		Transform inv = t->inverse();
		inv.copy_matrix_into_array(m);
		image->bind_cuda_texture();
		EMDataForCuda* tmp = emdata_transform_cuda(m,image->get_xsize(),image->get_ysize(),image->get_zsize());
		image->unbind_cuda_texture();
		delete [] m;
		if (tmp == 0) throw;
		image->set_gpu_rw_data(tmp->data,tmp->nx,tmp->ny,tmp->nz);
		free(tmp);
	}
#endif
	if ( use_cpu ) {
		float* des_data = transform(image,*t);
		image->set_data(des_data,image->get_xsize(),image->get_ysize(),image->get_zsize());
		image->update();
	}
	float scale = t->get_scale();
	if (scale != 1.0) {
		image->scale_pixel(1.0f/scale);
//		update_emdata_attributes(image,image->get_attr_dict(),scale);
	}

	if(t) {delete t; t=0;}

	EXITFUNC;
}


void IntTranslateProcessor::assert_valid_aspect(const vector<int>& translation, const EMData* const image) const {
	if (translation.size() == 0 ) throw InvalidParameterException("You must specify the trans argument");
}

Region IntTranslateProcessor::get_clip_region(vector<int>& translation, const EMData* const image) const {
	unsigned int dim = static_cast<unsigned int> (image->get_ndim());

	if ( translation.size() != dim ) {
		for(unsigned int i = translation.size(); i < dim; ++i ) translation.push_back(0);
	}

	Region clip_region;
	if (dim == 1) {
		clip_region = Region(-translation[0],image->get_xsize());
	} else if ( dim == 2 ) {
		clip_region = Region(-translation[0],-translation[1],image->get_xsize(),image->get_ysize());
	} else if ( dim == 3 ) {
		clip_region = Region(-translation[0],-translation[1],-translation[2],image->get_xsize(),image->get_ysize(),image->get_zsize());
	} else throw ImageDimensionException("Only 1,2 and 3D images are supported");

	return clip_region;
}

void IntTranslateProcessor::process_inplace(EMData* image) {

	vector<int> translation = params.set_default("trans",vector<int>() );


	assert_valid_aspect(translation,image);

	Region clip_region = get_clip_region(translation,image);

	image->clip_inplace(clip_region,0);
	// clip_inplace does the update!
}

EMData* IntTranslateProcessor::process(const EMData* const image) {

	vector<int> translation = params.set_default("trans",vector<int>() );

	assert_valid_aspect(translation,image);

	Region clip_region = get_clip_region(translation,image);

	return image->get_clip(clip_region,0);
	// clip_inplace does the update!
}


void ScaleTransformProcessor::process_inplace(EMData* image) {
	int ndim = image->get_ndim();
	if (ndim != 2 && ndim != 3) throw UnexpectedBehaviorException("The Scale Transform processors only works for 2D and 3D images");

	if ( image->get_xsize() != image->get_ysize()) {
		throw ImageDimensionException("x size and y size of image do not match. This processor only works for uniformly sized data");
	}
	if ( ndim == 3 ) {
		if ( image->get_xsize() != image->get_zsize()) {
		throw ImageDimensionException("x size and z size of image do not match. This processor only works for uniformly sized data");
		}
	}

	float scale = params.set_default("scale",0.0f);
	if (scale <= 0.0f) throw InvalidParameterException("The scale parameter must be greater than 0");

	int clip = 0;

	if(params.has_key("clip"))
	{
		clip = params["clip"];
		if (clip < 0) throw InvalidParameterException("The clip parameter must be greater than 0"); // If it's zero it's not used
	}
	else
	{
		clip = (int)(scale*image->get_xsize());
	}

	Region r;
	if (ndim == 3) {
		 r = Region( (image->get_xsize()-clip)/2, (image->get_xsize()-clip)/2, (image->get_xsize()-clip)/2,clip, clip,clip);
	} else { // ndim == 2 guaranteed by check at beginning of function
		 r = Region( (image->get_xsize()-clip)/2, (image->get_xsize()-clip)/2, clip, clip);
	}

	if (scale > 1) {
		if ( clip != 0) {
			image->clip_inplace(r);
		}
		Transform t;
		t.set_scale(scale);
		image->process_inplace("xform",Dict("transform",&t));
	} else if (scale < 1) {
		Transform t;
		t.set_scale(scale);
		image->process_inplace("xform",Dict("transform",&t));
		if ( clip != 0) {
			image->clip_inplace(r);
		}
	} else {
		if ( clip != 0) {
			image->clip_inplace(r);
		}
	}
}

EMData* ScaleTransformProcessor::process(const EMData* const image) {
	int ndim = image->get_ndim();
	if (ndim != 2 && ndim != 3) throw UnexpectedBehaviorException("The Scale Transform processors only works for 2D and 3D images");

	if ( image->get_xsize() != image->get_ysize()) {
		throw ImageDimensionException("x size and y size of image do not match. This processor only works for uniformly sized data");
	}
	if ( ndim == 3 ) {
		if ( image->get_xsize() != image->get_zsize()) {
		throw ImageDimensionException("x size and z size of image do not match. This processor only works for uniformly sized data");
		}
	}

	float scale = params.set_default("scale",0.0f);
	if (scale <= 0.0f) throw InvalidParameterException("The scale parameter must be greater than 0");

	int clip = 0;

	if(params.has_key("clip"))
	{
		clip = params["clip"];
		if (clip < 0) throw InvalidParameterException("The clip parameter must be greater than 0"); // If it's zero it's not used
	}
	else
	{
		clip = (int)(scale*image->get_xsize());
	}

	Region r;
	if (ndim == 3) {
		 r = Region( (image->get_xsize()-clip)/2, (image->get_xsize()-clip)/2, (image->get_xsize()-clip)/2,clip, clip,clip);
	} else { // ndim == 2 guaranteed by check at beginning of function
		 r = Region( (image->get_xsize()-clip)/2, (image->get_xsize()-clip)/2, clip, clip);
	}

	EMData* ret = 0;
	if (scale > 1) {
		if ( clip != 0) {
			ret = image->get_clip(r);
		}
		Transform t;
		t.set_scale(scale);
		if (ret != 0) {
			ret->process_inplace("xform",Dict("transform",&t));
		} else {
			ret = image->process("xform",Dict("transform",&t));
		}
	} else if (scale < 1) {
		Transform t;
		t.set_scale(scale);
		ret = image->process("xform",Dict("transform",&t));
		if ( clip != 0) {
			ret->clip_inplace(r);
		}
	} else {
		if ( clip != 0) {
			ret = image->get_clip(r);
		} else {
			ret = image->copy();
		}
	}
	return ret;

}

void Rotate180Processor::process_inplace(EMData* image) {
	ENTERFUNC;


	if (image->get_ndim() != 2) {
		throw ImageDimensionException("2D only");
	}

#ifdef EMAN2_USING_CUDA
	if (image->gpu_operation_preferred() ) {
		EMDataForCuda tmp = image->get_data_struct_for_cuda();
		emdata_rotate_180(&tmp);
		image->gpu_update();
		EXITFUNC;
		return;
	}
#endif

	float *d = image->get_data();
	int nx = image->get_xsize();
	int ny = image->get_ysize();

	// x and y offsets are used to handle even vs odd cases
	int x_offset = 0;
	if (nx % 2 == 1) x_offset=1;
	int y_offset = 0;
	if (ny % 2 == 1) y_offset=1;

	bool stop = false;
	for (int x = 1; x <= (nx/2+x_offset); x++) {
		int y = 0;
		for (y = 1; y < (ny+y_offset); y++) {
			if (x == (nx / 2+x_offset) && y == (ny / 2+y_offset)) {
				stop = true;
				break;
			}
			int i = (x-x_offset) + (y-y_offset) * nx;
			int k = nx - x + (ny - y) * nx;

			float t = d[i];
			d[i] = d[k];
			d[k] = t;
		}
		if (stop) break;
	}

	/* Here we guard against irregularites that occur at the boundaries
	 * of even dimensioned images. The basic policy is to replace the pixel
	 * in row 0 and/or column 0 with those in row 1 and/or column 1, respectively.
	 * The pixel at 0,0 is replaced with the pixel at 1,1 if both image dimensions
	 * are even. FIXME - it may be better to use an average at the corner, in
	 * this latter case, using pixels (1,1), (0,1) and (1,0). I am not sure. (dsawoolford)
	*/
	if (x_offset == 0) {
		for (int y = 0; y < ny; y++) {
			image->set_value_at_fast(0,y,image->get_value_at(1,y));
		}
	}

	if (y_offset == 0) {
		for (int x = 0; x < nx; x++) {
			image->set_value_at_fast(x,0,image->get_value_at(x,1));
		}
	}

	if (y_offset == 0 && x_offset == 0) {
		image->set_value_at_fast(0,0,image->get_value_at(1,1));
	}

	image->update();
	EXITFUNC;
}


void ClampingProcessor::process_inplace( EMData* image )
{

	if ( image->is_complex() ) throw ImageFormatException("Error: clamping processor does not work on complex images");

	float min = params.set_default("minval",default_min);
	float max = params.set_default("maxval",default_max);
	bool tomean = params.set_default("tomean",false);
	float new_min_vals = min;
	float new_max_vals = max;
	if (tomean){

		new_min_vals = image->get_attr("mean");
		new_max_vals = image->get_attr("mean");
	}

	// Okay, throwing such an error is probably overkill - but atleast the user will get a loud message
	// saying what went wrong.
	if ( max < min ) throw InvalidParameterException("Error: minval was greater than maxval, aborting");

	size_t size = image->get_size();
	for(size_t i = 0; i < size; ++i )
	{
		float * data = &image->get_data()[i];
		if ( *data < min ) *data = new_min_vals;
		else if ( *data > max ) *data = new_max_vals;
	}
	image->update();
}

void TestTomoImage::insert_rectangle( EMData* image, const Region& region, const float& value, const Transform& t3d )
{
	int startx = (int)region.origin[0] - (int)region.size[0]/2;
	int starty = (int)region.origin[1] - (int)region.size[1]/2;
	int startz = (int)region.origin[2] - (int)region.size[2]/2;

	int endx  = (int)region.origin[0] + (int)region.size[0]/2;
	int endy  = (int)region.origin[1] + (int)region.size[1]/2;
	int endz  = (int)region.origin[2] + (int)region.size[2]/2;

	if ( ! t3d.is_identity() ) {
		float xt, yt, zt;
		for ( float z = (float)startz; z < (float)endz; z += 0.25f ) {
			for ( float y = (float)starty; y < (float)endy; y += 0.25f ) {
				for ( float x = (float)startx; x < (float)endx; x += 0.25f ) {
					xt = (float) x - region.origin[0];
					yt = (float) y - region.origin[1];
					zt = (float) z - region.origin[2];
					Vec3f v((float)xt,(float)yt,(float)zt);
					v = t3d*v;
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

void TestTomoImage::process_inplace( EMData* image )
{
	//float nx = 240;
	//float ny = 240;
	//float nz = 60;

	//image->set_size((int)nx,(int)ny,(int)nz);
	float nx = (float) image->get_xsize();
	float ny = (float) image->get_ysize();
	float nz = (float) image->get_zsize();

	// This increment is used to simplified positioning
	// It's an incremental factor that matches the grid size of the paper
	// that I drew this design on before implementing it in code
	float inc = 1.0f/22.0f;
	float xinc = inc;
	float yinc = inc;
	float zinc = inc;

	Dict d;
	d["a"] = (float) .4*nx+3;
	d["b"] = (float) .4*ny+3;
	d["c"] = (float) .4*nz+3;
	d["fill"] = 0.2;
	image->process_inplace("testimage.ellipsoid",d);

	d["a"] = (float) .4*nx;
	d["b"] = (float) .4*ny;
	d["c"] = (float) .4*nz;
	d["fill"] = 0.1;
	image->process_inplace("testimage.ellipsoid",d);

	// Center x, center z, bottom y ellipsoids that grow progessively smaller
	{
		Transform t;
		t.set_trans(0.,ny*4.0f*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float) 2.*xinc*nx;
		d["b"] = (float)0.5*yinc*ny;
		d["c"] = (float) 1.*zinc*nz;
		d["fill"] = 0.3;
		image->process_inplace("testimage.ellipsoid",d);
	}

	{
		Transform t;
		t.set_trans(0.,ny*5.5f*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float) 1.5*xinc*nx;
		d["b"] = (float)0.5*yinc*ny;
		d["c"] = (float) 1.*zinc*nz;
		d["fill"] = 0.0;
		image->process_inplace("testimage.ellipsoid",d);
	}
	{
		Transform t;
		t.set_trans(0.,ny*7*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float) 1.*xinc*nx;
		d["b"] = (float)0.5*yinc*ny;
		d["c"] = (float) 1.*zinc*nz;
		d["fill"] = 0.3;
		image->process_inplace("testimage.ellipsoid",d);
	}


	{
		Transform t;
		t.set_trans(0.,ny*8.5f*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float) .75*xinc*nx;
		d["b"] = (float)0.5*yinc*ny;
		d["c"] = (float) 1.*zinc*nz;
		d["fill"] = 0.0;
		image->process_inplace("testimage.ellipsoid",d);
	}

	// Center x, center z, bottom y ellipsoids that grow progessively smaller
	{
		Transform t;
		t.set_trans(0.,ny*18*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float) 2*xinc*nx;
		d["b"] = (float)0.5*yinc*ny;
		d["c"] = (float) 1.*zinc*nz;
		d["fill"] = 0.3;
		image->process_inplace("testimage.ellipsoid",d);
	}

	{
		Transform t;
		t.set_trans(0.,ny*16.5f*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float) 1.5*xinc*nx;
		d["b"] = (float)0.5*yinc*ny;
		d["c"] = (float) 1.*zinc*nz;
		d["fill"] = 0.3;
		image->process_inplace("testimage.ellipsoid",d);
	}

	{
		Transform t;
		t.set_trans(0.,ny*15*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float) 1*xinc*nx;
		d["b"] = (float)0.5*yinc*ny;
		d["c"] = (float) 1.*zinc*nz;
		d["fill"] = 0.3f;
		image->process_inplace("testimage.ellipsoid",d);
	}

	{
		Transform t;
		t.set_trans(0.,ny*13.5f*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float).75*xinc*nx;
		d["b"] = (float)0.5*yinc*ny;
		d["c"] = (float) 1.*zinc*nz;
		d["fill"] = 0.3;
		image->process_inplace("testimage.ellipsoid",d);
	}

	// Left ellipsoids from the bottom up
	{

		Transform t;
		t.set_trans(nx*6*xinc-nx/2,ny*5*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float)1*xinc*nx;
		d["b"] = (float).75*yinc*ny;
		d["c"] = (float) .75*zinc*nz;
		d["fill"] = 0.25;
		image->process_inplace("testimage.ellipsoid",d);
	}

	{
		Transform t;
		t.set_trans(nx*6*xinc-nx/2,ny*7*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float)1.5*xinc*nx;
		d["b"] = (float).75*yinc*ny;
		d["c"] = (float) .75*zinc*nz;
		d["fill"] = 0.25;
		image->process_inplace("testimage.ellipsoid",d);
	}

	{
		Transform t;
		t.set_trans(nx*6*xinc-nx/2,ny*9*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float)2*xinc*nx;
		d["b"] = (float).75*yinc*ny;
		d["c"] = (float) .75*zinc*nz;
		d["fill"] = 0.25;
		image->process_inplace("testimage.ellipsoid",d);
	}

	{
		Transform t;
		t.set_trans(nx*6*xinc-nx/2,ny*11*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float)2.5*xinc*nx;
		d["b"] = (float).75*yinc*ny;
		d["c"] = (float) 1*zinc*nz;
		d["fill"] = 0.25;
		image->process_inplace("testimage.ellipsoid",d);
	}

	{
		Transform t;
		t.set_trans(nx*6*xinc-nx/2,ny*13*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float) 3*xinc*nx;
		d["b"] = (float).75*yinc*ny;
		d["c"] = (float) 1*zinc*nz;
		d["fill"] = 0.25;
		image->process_inplace("testimage.ellipsoid",d);
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
		Transform t3d(Dict("type","eman","az",(float)-25.0));
		insert_rectangle(image, region, 0.4f, t3d);
	}

	// Rotated ellipsoids
	{
		Transform t;
		t.set_trans(nx*6.8f*xinc-nx/2,ny*16*yinc-ny/2,0);
		Dict rot;
		rot["type"] = "eman";
		rot["az"] = 43.0f;
		t.set_rotation(rot);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float) 1.5*xinc*nx;
		d["b"] = (float) .5*yinc*ny;
		d["c"] = (float) .5*zinc*nz;
		d["fill"] = 0.2;
		image->process_inplace("testimage.ellipsoid",d);
	}
	{
		Transform t;
		t.set_trans(nx*7.2f*xinc-nx/2,ny*16*yinc-ny/2,0);
		Dict rot;
		rot["type"] = "eman";
		rot["az"] = 135.0f;
		t.set_rotation(rot);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float) 1.5*xinc*nx;
		d["b"] = (float) .5*yinc*ny;
		d["c"] = (float) .5*zinc*nz;
		d["fill"] = 0.3;
		image->process_inplace("testimage.ellipsoid",d);
	}

	// Dense small ellipsoids
	{
		Transform t;
		t.set_trans(nx*3.5f*xinc-nx/2,ny*8*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float) .5*xinc*nx;
		d["b"] = (float) .5*yinc*ny;
		d["c"] = (float) .5*zinc*nz;
		d["fill"] = 2.05;
		image->process_inplace("testimage.ellipsoid",d);

		t.set_trans(nx*8*xinc-nx/2,ny*18*yinc-ny/2,0);
		image->process_inplace("testimage.ellipsoid",d);

		t.set_trans(nx*14*xinc-nx/2,ny*18.2f*yinc-ny/2,0);
		image->process_inplace("testimage.ellipsoid",d);

		t.set_trans(nx*18*xinc-nx/2,ny*14*yinc-ny/2,0);
		image->process_inplace("testimage.ellipsoid",d);

		t.set_trans(nx*17*xinc-nx/2,ny*7.5f*yinc-ny/2,0);
		image->process_inplace("testimage.ellipsoid",d);
	}


	// Dense small rectangles
	{
		Region region(nx*18.*inc,ny*11.5*inc,nz/2.,1.*inc*nx,1.*inc*ny,1.*inc*nz);
		Transform t3d(Dict("type","eman","az",(float)45.0));
		insert_rectangle(image, region, 1.45f, t3d);
	}
	{
		Region region(nx*3.*inc,ny*10.5*inc,nz/2.,1.*inc*nx,1.*inc*ny,1.*inc*nz);
		Transform t3d(Dict("type","eman","az",(float)45.0));
		insert_rectangle(image, region, 1.45f, t3d);
	}

	// Insert small cluster of spheres
	{
		Transform t;
		t.set_trans(nx*14*xinc-nx/2,ny*7.5f*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float) .5*xinc*nx;
		d["b"] = (float) .5*yinc*ny;
		d["c"] = (float) .5*zinc*nz;
		d["fill"] = .35;
		image->process_inplace("testimage.ellipsoid",d);
	}
	{
		Transform t;
		t.set_trans(nx*15*xinc-nx/2,ny*7.5f*yinc-ny/2,0);
		Dict d;
		d["transform"] = &t;
		d["a"] = (float) .25*xinc*nx;
		d["b"] = (float) .25*yinc*ny;
		d["c"] = (float) .25*zinc*nz;
		d["fill"] = .35;
		image->process_inplace("testimage.ellipsoid",d);

		t.set_trans(nx*13.5f*xinc-nx/2,ny*6.5f*yinc-ny/2,0);
		image->process_inplace("testimage.ellipsoid",d);

		t.set_trans(nx*14.5f*xinc-nx/2,ny*6.5f*yinc-ny/2,0);
		image->process_inplace("testimage.ellipsoid",d);

		t.set_trans(nx*15.5f*xinc-nx/2,ny*6.5f*yinc-ny/2,0);
		image->process_inplace("testimage.ellipsoid",d);

		t.set_trans(nx*14*xinc-nx/2,ny*5.5f*yinc-ny/2,0);
		image->process_inplace("testimage.ellipsoid",d);

		t.set_trans(nx*14*xinc-nx/2,ny*5.5f*yinc-ny/2,0);
		image->process_inplace("testimage.ellipsoid",d);

		t.set_trans(nx*15*xinc-nx/2,ny*5.5f*yinc-ny/2,0);
		image->process_inplace("testimage.ellipsoid",d);

		t.set_trans(nx*16*xinc-nx/2,ny*5.5f*yinc-ny/2,0);
		image->process_inplace("testimage.ellipsoid",d);

		t.set_trans(nx*14.5f*xinc-nx/2,ny*4.5f*yinc-ny/2,0);
		image->process_inplace("testimage.ellipsoid",d);

		t.set_trans(nx*15.5f*xinc-nx/2,ny*4.5f*yinc-ny/2,0);
		image->process_inplace("testimage.ellipsoid",d);
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
	float nbins = (float)params.set_default("nbins",default_bins);
	bool debug = params.set_default("debug",false);

	vector<int> debugscores;
	if ( debug ) {
		debugscores = vector<int>((int)nbins, 0);
	}

	if ( nbins < 0 ) throw InvalidParameterException("nbins must be greater than 0");

	float bin_width = (max-min)/nbins;
	float bin_val_offset = bin_width/2.0f;

	size_t size = image->get_size();
	float* dat = image->get_data();

	for(size_t i = 0; i < size; ++i ) {
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
	//bool complexflag = false;
	EMData* null = 0;
	EMData* with = params.set_default("with", null);
	if ( with == NULL ) throw InvalidParameterException("Error - the image required for the convolution is null");

	EMData* newimage = fourierproduct(image, with, CIRCULANT, CONVOLUTION, false);

	float* orig = image->get_data();
	float* work = newimage->get_data();
	int nx  = image->get_xsize();
	int ny  = image->get_ysize();
	int nz  = image->get_zsize();
	memcpy(orig,work,nx*ny*nz*sizeof(float));
	image->update();

	delete newimage;
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
	image->process_inplace("math.convolution", conv_parms);
	image->process_inplace("xform.phaseorigin.tocenter");

	delete e;
}


void TomoTiltAngleWeightProcessor::process_inplace( EMData* image )
{
	bool fim = params.set_default("angle_fim", false);
	float alt;
	if ( fim ) {
		alt = image->get_attr("euler_alt");
	}
	else alt = params.set_default("angle", 0.0f);

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
		Transform* t = (Transform*)image->get_attr("xform.projection");
		Dict d = t->get_params("eman");
		alt = (float) d["alt"];
		if(t) {delete t; t=0;}
	}
	else alt = params.set_default("angle", 0.0f);


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
	image->process_inplace("math.convolution", conv_parms);

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
	image->process_inplace("math.convolution", conv_parms);

	delete e;
}


#ifdef EMAN2_USING_CUDA

/* CLASS MPICUDA_kmeans processor
 *
 */
#include "sparx/cuda/cuda_mpi_kmeans.h"
MPICUDA_kmeans::MPICUDA_kmeans() {
    h_IM = NULL;
    h_AVE = NULL;
    h_ASG = NULL;
    h_dist = NULL;
    h_AVE2 = NULL;
    h_im2 = NULL;
    h_NC = NULL;
    params = NULL;
    d_im = NULL;
    d_AVE = NULL;
    d_dist = NULL;
}

MPICUDA_kmeans::~MPICUDA_kmeans() {
    if (h_IM) delete h_IM;
    if (h_ASG) delete h_ASG;
    if (h_AVE) delete h_AVE;
    if (h_dist) delete h_dist;
    if (h_AVE2) delete h_AVE2;
    if (h_im2) delete h_im2;
    if (h_NC) delete h_NC;
    if (params) delete params;
}

#include "sparx/cuda/cuda_mpi_kmeans.h"
int MPICUDA_kmeans::setup(int extm, int extN, int extn, int extK, int extn_start) {
    m = extm;				// number of pixels per image
    N = extN;				// Total number of images
    n = extn;                           // Number of images used locally
    K = extK;				// number of classes
    n_start = extn_start;               // Starting point to local images
    size_IM = m * N;                    // nb elements in IM
    size_im = m * n;                    // nb elements in im
    size_AVE = m * K;                   // nb elements in ave
    size_dist = n * K;                  // nb elements in dist
    ite = 0;                            // init nb of iterations
    BLOCK_SIZE = 512;
    NB = size_dist / BLOCK_SIZE;
    ins_BLOCK = NB * BLOCK_SIZE;
    // Host memory allocation for images
    h_IM = (float*)malloc(size_IM * sizeof(float));
    if (h_IM == 0) return 1;
    h_im = &h_IM[n_start * m]; // for local images
    // Host memory allocation for the averages
    h_AVE = (float*)malloc(size_AVE * sizeof(float));
    if (h_AVE == 0) return 1;
    // Host memory allocation for all assignment
    h_ASG = (unsigned short int*)malloc(N * sizeof(unsigned short int));
    if (h_ASG == 0) return 1;
    h_asg = &h_ASG[n_start]; // for local assignment
    // allocate the memory for the sum squared of averages
    h_AVE2 = (float*)malloc(K * sizeof(float));
    if (h_AVE2 == 0) return 1;
    // allocate the memory for the sum squared of images
    h_im2 = (float*)malloc(n * sizeof(float));
    if (h_im2 == 0) return 1;
    // allocate the memory for the distances
    h_dist = (float*)malloc(size_dist * sizeof(float));
    if (h_dist == 0) return 1;
    // allocate the memory for the number of images per class
    h_NC = (unsigned int*)malloc(K * sizeof(unsigned int));
    if (h_NC == 0) return 1;
    // allocate the memory to parameters
    params = (int*)malloc(8 * sizeof(int));
    if (params == 0) return 1;
    params[0] = n;
    params[1] = m;
    params[2] = K;
    params[3] = 0;          // Reserve for flag_stop_iteration
    params[4] = 0;          // Reserve for ct_im_mv (debug)
    params[5] = BLOCK_SIZE; // Size of threads block (512)
    params[6] = NB;         // Number of blocks which fit with BLOCK_SIZE
    params[7] = ins_BLOCK;  // Number of blocks remaining

    return 0;
}

// add image pre-process by Util.compress_image_mask
void MPICUDA_kmeans::append_flat_image(EMData* im, int pos) {
    for (int i = 0; i < m ; ++i) h_IM[pos * m + i] = (*im)(i);
}

// cuda init mem device, cublas (get device ptr)
int MPICUDA_kmeans::init_mem(int numdev) {
    int stat = 1;
    float** hd_AVE = NULL;
    float** hd_im = NULL;
    float** hd_dist = NULL;
    hd_AVE = &d_AVE;
    hd_im = &d_im;
    hd_dist = &d_dist;
    stat = cuda_mpi_init(h_im, hd_im, hd_AVE, hd_dist, size_im, size_AVE, size_dist, numdev);
    //printf("C++ get this pointer for d_AVE: %p\n", d_AVE);
    //printf("C++ get this pointer for d_im: %p\n", d_im);
    //printf("C++ get this pointer for d_dist: %p\n", d_dist);
    return stat;
}

// precalculate im2
void MPICUDA_kmeans::compute_im2() {
    for (int i = 0; i < n; i++) {
	h_im2[i] = 0.0f;
	for (int j = 0; j < m; j++) h_im2[i] += (h_im[i * m + j] * h_im[i * m + j]);
    }
}

// init randomly the first assignment
int MPICUDA_kmeans::random_ASG(long int rnd) {
    srand(rnd);
    int ret = 20;
    int flag = 0;
    int i, k;
    for (k = 0; k < K; k++) h_NC[k] = 0;
    while (ret > 0) {
	ret--;
	for (i = 0; i < N; i++) {
	    h_ASG[i] = rand() % K;
	    h_NC[h_ASG[i]]++;
	}
	flag = 1;
	k = K;
	while (k > 0 && flag) {
	    k--;
	    if (h_NC[k] <= 1) {
		flag = 0;
		if (ret == 0) {
		    //printf("Erreur randomize assignment\n");
		    return -1;
		}
		for (k = 0; k < K; k++) h_NC[k] = 0;
	    }
	if (flag == 1) ret = 0;
	}
    }
    return 0;
}

// get back the assignment
vector <int> MPICUDA_kmeans::get_ASG() {
    vector <int> ASG(h_ASG, &h_ASG[N]);
    return ASG;
}

// get back the global assignment
vector <int> MPICUDA_kmeans::get_asg() {
    vector <int> asg(h_asg, &h_asg[n]);
    return asg;
}

// compute NC from ASG
void MPICUDA_kmeans::compute_NC() {
    for (int i = 0; i < K; ++i) h_NC[i] = 0;
    for (int i = 0; i < N; ++i) h_NC[h_ASG[i]]++;
}

// get NC
vector <int> MPICUDA_kmeans::get_NC() {
    vector <int> NC(h_NC, &h_NC[K]);
    return NC;
}

// set a new global assignment
void MPICUDA_kmeans::set_ASG(const vector <int>& ASG) {
    for (int i = 0; i < N ; ++i) h_ASG[i] = ASG[i];
}

// set number of objects per group
void MPICUDA_kmeans::set_NC(const vector <int>& NC) {
    for (int i = 0; i < K; ++i) h_NC[i] = NC[i];
}

// get back some information (ite and T0)
int MPICUDA_kmeans::get_ct_im_mv() {
    return params[4]; // ct_im_mov
}

// set the value of T
void MPICUDA_kmeans::set_T(float extT) {
    T = extT;
}

// get the T value
float MPICUDA_kmeans::get_T() {
    return T;
}

// compute ave and ave2
void MPICUDA_kmeans::compute_AVE() {
    float buf = 0.0f;
    int i, j, c, d, ind;
    // compute the averages according ASG
    for (i = 0; i < size_AVE; ++i) h_AVE[i] = 0.0f;                          // set averages to zero
    for (i = 0; i < N; ++i) {
	c = h_ASG[i] * m;
	d = i * m;
	for (j = 0; j < m; ++j) h_AVE[c + j] += h_IM[d + j];                 // accumulate images
    }
    for (i = 0; i < K; i++) {
	buf = 0.0f;
	for (j = 0 ; j < m; j++) {
	    ind = i * m + j;
	    h_AVE[ind] /= (float)h_NC[i];                                    // compute average
	    buf += (h_AVE[ind] * h_AVE[ind]);                                // do sum squared AVE
	}
	h_AVE2[i] = buf;
    }
}

// set averages
void MPICUDA_kmeans::set_AVE(EMData* im, int pos) {
    for (int i = 0; i < m ; ++i) h_AVE[pos * m + i] = (*im)(i);
}

// get back the averages
vector<EMData*> MPICUDA_kmeans::get_AVE() {
    vector<EMData*> ave(K);
    for (int k = 0; k < K; ++k) {
	EMData* im = new EMData();
	im->set_size(m, 1, 1);
	float *ptr = im->get_data();
	for (int i = 0; i < m; ++i) {ptr[i] =  h_AVE[k * m + i];}
	ave[k] = im->copy();
	delete im;
    }
    return ave;
}

// k-means one iteration
int MPICUDA_kmeans::one_iter() {
    int status = cuda_mpi_kmeans(h_AVE, d_AVE, h_dist, d_dist, d_im, h_im2, h_AVE2, h_asg, h_NC, params);
    ite++;
    return status;
}

// k-means SA one iteration
int MPICUDA_kmeans::one_iter_SA() {
    int status = cuda_mpi_kmeans_SA(h_AVE, d_AVE, h_dist, d_dist, d_im, h_im2, h_AVE2, h_asg, h_NC, T, params);
    ite++;
    return status;
}

// compute ji
vector <float> MPICUDA_kmeans::compute_ji() {
    int status = cuda_mpi_dist(h_AVE, d_AVE, h_dist, d_dist, d_im, n, K, m);
    vector <float> ji(K);
    int i;
    if (status != 0) {
	for (i = 0; i < K; ++i) ji[i] = -1.0;
	return ji;
    }
    for (i = 0; i < n; ++i) ji[h_asg[i]] += (h_im2[i] + h_AVE2[h_asg[i]] - 2 * h_dist[i * K + h_asg[i]]);
    return ji;
}

//
vector <float> MPICUDA_kmeans::compute_criterion(const vector <float>& Ji) {
    float buf = 0.0f;
    float Je = 0.0f;
    float Tr_AVE = 0.0f;
    float v_max = 0.0f;
    float* S_AVE2 = (float*)calloc(m, sizeof(float));
    float* S_AVE = (float*)calloc(m, sizeof(float));
    vector <float> crit(4);
    int i, j, k;
    // Je
    for (i = 0; i < K; ++i) Je += (Ji[i] / float(m));
    crit[0] = Je;
    // trace ave
    for (i = 0; i < K; ++i) {
	for (j = 0; j < m; ++j) {
	    S_AVE[j] += h_AVE[i * m + j];
	    S_AVE2[j] += (h_AVE[i * m + j] * h_AVE[i * m +j]);
	}
    }
    buf = 1 / (float)K;
    for (i = 0; i < m; ++i) Tr_AVE += (buf * (S_AVE2[i] - buf * S_AVE[i] * S_AVE[i]));
    // Coleman
    crit[1] = Tr_AVE * Je;
    // Harabasz
    crit[2] = (Tr_AVE * (float)(N - K)) / (Je * (float)(K - 1));
    // Davies-Bouldin
    for (i = 0; i < K; ++i) {
	v_max = 0.0f;
	for (j = 0; j < K; ++j) {
	    if (j != i) {
		buf = 0.0f;
		for (k = 0; k < m; ++k) buf += ((h_AVE[j * m + k] - h_AVE[i * m + k]) * (h_AVE[j * m + k] - h_AVE[i * m + k]));
		buf = (Ji[i] / (float)h_NC[i] + Ji[j] / (float)h_NC[j]) * ((float)m / buf);
	    }
	    if (buf > v_max) v_max = buf;
	}
	crit[3] += v_max;
    }
    crit[3] /= (float)K;
    free(S_AVE);
    free(S_AVE2);
    return crit;
}

// shutdown cublas and release device mem
int MPICUDA_kmeans::shutdown() {
    return cuda_mpi_shutdown(d_im, d_AVE, d_dist);
}
//// END OF CUDA KMEANS /////////////////////////

void CudaMultProcessor::process_inplace(EMData* image) {
	float val = params.set_default("scale",(float) 1.0);
	EMDataForCuda tmp = image->get_data_struct_for_cuda();
	emdata_processor_mult(&tmp,val);
	image->gpu_update();
}


void CudaCorrelationProcessor::process_inplace(EMData* image) {
	EMData* with = params.set_default("with",(EMData*)0);
	if (with == 0) throw InvalidParameterException("You must supply the with parameter, and it must be valid. It is NULL.");

	EMDataForCuda left = image->get_data_struct_for_cuda();
	with->bind_cuda_texture(false);
	emdata_processor_correlation_texture(&left,1);
	image->gpu_update();

}

#endif //EMAN2_USING_CUDA

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

float ModelHelixProcessor::radprofile(float r, int type)
//Ross Coleman: modified from EMAN1 Cylinder.C by Wen Jiang
{
	// r in angstrom
	double ret = 0;//iterations ==> rounding is problematic with float types if 15 < r < 20 and really bad if 20 < r < 30
	if (type == 0) { // pure Gaussian falloff
		r /= 2;
		ret = exp(-r * r);
	} else if (type == 1) { // pure Gaussian falloff + negative dip, so mean is 0
		r /= 2;
		ret = (1 - r * r / 4) * exp(-r * r / 4);
	} else if (type == 2) {
		// polynomial fitting to the radial profile of real helix density
		// f=a0*x^n+a1+x^(n-1)+ ... +a[n-1]*x+an
		//float an[11]={2.847024584977009e-10,-3.063997224364090e-08,1.418801040660860e-06,-3.678676414383996e-05,5.804871622801710e-04,-5.640340018430164e-03,3.208802421493864e-02,-9.068475313823952e-02,7.097329559749284e-02,-9.993347339658298e-02,1.000000000000000e+00};

		// now the fitting to the original profile
		if (r >= 12.2)
			return 0; //We don't want that part of the polynomial --> goes way below zero
		static float an[15] = { -3.9185246832229140e-16f,
				3.3957205298900993e-14f, 2.0343351971222658e-12f,
				-4.4935965816879751e-10f, 3.0668169835080933e-08f,
				-1.1904544689091790e-06f, 2.9753088549414953e-05f,
				-4.9802112876220150e-04f, 5.5900917825309360e-03f,
				-4.0823714462925299e-02f, 1.8021733669148599e-01f,
				-4.0992557296268717e-01f, 3.3980328566901458e-01f,
				-3.6062024812411908e-01f, 1.0000000000000000e+00f };

		ret = an[0];
		for (int i = 1; i < 15; i++) {
			ret = ret * r + an[i];
		}
	}
	return (float)ret;
}

void ModelEMCylinderProcessor::process_inplace(EMData * in)
//Ross Coleman: modified from EMAN1 Cylinder.C by Wen Jiang
{
	// synthesize model alpha helix, len is Angstrom, default to 2 turns
	//The helical axis is parallel to the z axis.
	EMData * cyl = in;
	int nx = cyl->get_xsize();
	int ny = cyl->get_ysize();
	int nz = cyl->get_zsize();

	int type = params.set_default("type", 2);
	float len = params.set_default("length", 16.2); //in angstroms
	int x0 = params.set_default("x0", -1); //in voxels -- default value changed a few lines down
	int y0 = params.set_default("y0", -1); //in voxels
	int z0 = params.set_default("z0", -1); //in voxels
	//TODO: check with Matt about default values

	if (x0 < 0 || x0 >= nx)
		x0 = nx / 2;
	if (y0 < 0 || y0 >= ny)
		y0 = ny / 2;
	if (z0 < 0 || z0 >= nz)
		z0 = nz / 2;

	float apix_x = cyl->get_attr("apix_x"); //TODO: Ask Matt if I correctly handled cases where apix_x != apix_y or apix_x != apix_z are not equal
	float apix_y = cyl->get_attr("apix_y");
	float apix_z = cyl->get_attr("apix_z");

	float * dat = cyl->get_data();
	int cyl_voxel_len = (int) (len / apix_z);
	int cyl_k_min = z0 - cyl_voxel_len / 2;
	int cyl_k_max = z0 + cyl_voxel_len / 2;

	int x, y;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++, dat++) {
				x = i - x0;//coordinate sys centered on cylinder
				y = j - y0;//coordinate sys centered on cylinder
				float radius = (float)hypot(x * apix_x, y * apix_y);
				if ((k > cyl_k_min) && (k < cyl_k_max))
					*dat += radprofile(radius, type); //pointer arithmetic for array done in loop
				//else
					//continue;

			}
		}
	}
}

void ApplyPolynomialProfileToHelix::process_inplace(EMData * in)
{
	EMData * cyl = in;
	int nx = cyl->get_xsize();
	int ny = cyl->get_ysize();
	int nz = cyl->get_zsize();
	float apix_x = cyl->get_attr("apix_x"); //TODO: Ask Matt if I correctly handled cases where apix_x != apix_y or apix_x != apix_z are not equal
	float apix_y = cyl->get_attr("apix_y");
	float apix_z = cyl->get_attr("apix_z");
	float lengthAngstroms = params["length"];
	int z0 = params.set_default("z0", -1); //in voxels

	if (z0 < 0 || z0 >= nz)
		z0 = nz / 2;

	int z_start = Util::round( z0 - 0.5*lengthAngstroms/apix_z );
	int z_stop = Util::round( z0 + 0.5*lengthAngstroms/apix_z );

	float * dat = cyl->get_data();
	double rho_x_sum, rho_y_sum, rho_sum, x_cm, y_cm, radius;

	for (int k = 0; k < nz; k++) //taking slices along z axis
	{
		rho_x_sum = rho_y_sum = rho_sum = 0; //Set to zero for a new slice

		if (k >= z_start && k <= z_stop)
		//Apply the radial profile only between z_start and z_stop on the z axis
		{
			//Calculating CM for the slice...
			for (int j = 0; j < ny; j++)
			{
				for (int i = 0; i < nx; i++, dat++)
				{
					rho_x_sum += (*dat)*i;
					rho_y_sum += (*dat)*j;
					rho_sum += *dat;
				}
			}
			if (rho_sum != 0)
			{
				dat -= nx*ny;//going back to the beginning of the dat array
				x_cm = rho_x_sum/rho_sum;
				y_cm = rho_y_sum/rho_sum;

				//Applying radial profile...
				for (int j=0; j<ny;j++)
				{
					for (int i=0;i<nx;i++,dat++)
					{
						radius = hypot( (i-x_cm)*apix_x, (j-y_cm)*apix_y );
						*dat = radprofile((float)radius, 2);//Type 2 is the polynomial radial profile.
					}
				}
			}
		}
		else
		//Clear the map, setting the density to zero everywhere else.
		{
			for (int j=0; j<ny; j++)
				for(int i=0; i<nx; i++)
				{
					*dat = 0;
					dat++;
				}
		}

	}
}

EMData* BinarySkeletonizerProcessor::process(EMData * image)
{
	using namespace wustl_mm::GraySkeletonCPP;
	using namespace wustl_mm::SkeletonMaker;

	Volume * vimage = new Volume(image);
	float threshold = params["threshold"];
	int min_curvew = params.set_default("min_curve_width", 4);
	int min_srfcw = params.set_default("min_surface_width", 4);
	bool mark_surfaces = params.set_default("mark_surfaces", true);
	Volume* vskel = VolumeSkeletonizer::PerformPureJuSkeletonization(vimage, "unused", static_cast<double>(threshold), min_curvew, min_srfcw);
	//VolumeSkeletonizer::CleanUpSkeleton(vskel, 4, 0.01f);
	if (mark_surfaces) {
		VolumeSkeletonizer::MarkSurfaces(vskel);
	}

	vskel->getVolumeData()->owns_emdata = false; //ensure the EMData object will remain when the Volume and its VolumeData object are freed
	EMData* skel = vskel->get_emdata();
	skel->update();
	return skel;
}

void BinarySkeletonizerProcessor::process_inplace(EMData * image)
{
	EMData* em_skel = this->process(image);
	//TODO: use memcpy and copy metadata explicitly
	*image = *em_skel; //Deep copy performed by EMData::operator=
	image->update();
	delete em_skel;
	em_skel = NULL;
	return;
}
