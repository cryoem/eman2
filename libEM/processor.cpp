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
#include "cmp.h"
#include "ctf.h"
#include "xydata.h"
#include "emdata.h"
#include "emassert.h"
#include "randnum.h"
#include "symmetry.h"
#include "averager.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_wavelet2d.h>
#include <gsl/gsl_multimin.h>
#include <algorithm>
#include <gsl/gsl_fit.h>
#include <ctime>

#ifdef __APPLE__
	typedef unsigned int uint;
#endif	//__APPLE__

#ifdef _WIN32
	typedef unsigned int uint;
#endif	//_WIN32

#ifdef EMAN2_USING_CUDA
//#include "cuda/cuda_util.h"
#include "cuda/cuda_processor.h"
#endif // EMAN2_USING_CUDA

using namespace EMAN;
using std::reverse;

const string SNREvalProcessor::NAME = "eval.maskedsnr";
const string AmpweightFourierProcessor::NAME = "filter.ampweight";
const string Axis0FourierProcessor::NAME = "filter.xyaxes0";
const string ConvolutionProcessor::NAME = "math.convolution";
const string XGradientProcessor::NAME = "math.edge.xgradient";
const string YGradientProcessor::NAME = "math.edge.ygradient";
const string ZGradientProcessor::NAME = "math.edge.zgradient";
const string Wiener2DAutoAreaProcessor::NAME = "filter.wiener2dauto";
const string Wiener2DFourierProcessor::NAME = "filter.wiener2d";
const string CtfSimProcessor::NAME = "math.simulatectf";
const string LinearRampFourierProcessor::NAME = "filter.linearfourier";
const string LoGFourierProcessor::NAME = "filter.LoG";
const string DoGFourierProcessor::NAME = "filter.DoG";
const string AzSharpProcessor::NAME = "filter.azimuthal.contrast";
const string HighpassAutoPeakProcessor::NAME = "filter.highpass.autopeak";
const string LinearRampProcessor::NAME = "eman1.filter.ramp";
const string AbsoluateValueProcessor::NAME = "math.absvalue";
const string FloorValueProcessor::NAME = "math.floor";
const string BooleanProcessor::NAME = "threshold.notzero";
const string KmeansSegmentProcessor::NAME = "segment.kmeans";
const string DistanceSegmentProcessor::NAME = "segment.distance";
const string WatershedProcessor::NAME = "segment.watershed";
const string RecipCarefullyProcessor::NAME = "math.reciprocal";
const string SubtractOptProcessor::NAME = "math.sub.optimal";
const string ValuePowProcessor::NAME = "math.pow";
const string ValueSquaredProcessor::NAME = "math.squared";
const string ValueSqrtProcessor::NAME = "math.sqrt";
const string ToZeroProcessor::NAME = "threshold.belowtozero";
const string AboveToZeroProcessor::NAME="threshold.abovetozero";
const string OutlierProcessor::NAME="threshold.outlier.localmean";
const string Rotate180Processor::NAME = "math.rotate.180";
const string TransformProcessor::NAME = "xform";
const string IntTranslateProcessor::NAME = "xform.translate.int";
const string ScaleTransformProcessor::NAME = "xform.scale";
const string ApplySymProcessor::NAME = "xform.applysym";
const string ClampingProcessor::NAME = "threshold.clampminmax";
const string NSigmaClampingProcessor::NAME = "threshold.clampminmax.nsigma";
const string ToMinvalProcessor::NAME = "threshold.belowtominval";
const string CutToZeroProcessor::NAME = "threshold.belowtozero_cut";
const string BinarizeProcessor::NAME = "threshold.binary";
//const string BinarizeAmpProcessor::NAME = "threshold.amp.binary";
const string BinarizeFourierProcessor::NAME = "threshold.binary.fourier";
const string CollapseProcessor::NAME = "threshold.compress";
const string LinearXformProcessor::NAME = "math.linear";
const string ExpProcessor::NAME = "math.exp";
const string FiniteProcessor::NAME = "math.finite";
const string RangeThresholdProcessor::NAME = "threshold.binaryrange";
const string SigmaProcessor::NAME = "math.sigma";
const string LogProcessor::NAME = "math.log";
const string MaskSharpProcessor::NAME = "mask.sharp";
const string MaskSoftProcessor::NAME = "mask.soft";
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
const string ZeroConstantProcessor::NAME = "mask.contract";	// This is broken, it never worked. Somebody didn't think it through properly
const string BoxMedianProcessor::NAME = "eman1.filter.median";
const string BoxSigmaProcessor::NAME = "math.localsigma";
const string NonConvexProcessor::NAME = "math.nonconvex";
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
const string BilateralProcessor::NAME = "filter.bilateral";
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
const string AddRandomNoiseProcessor::NAME = "math.addspectralnoise";
const string FourierToCornerProcessor::NAME = "xform.fourierorigin.tocorner";
const string FourierToCenterProcessor::NAME = "xform.fourierorigin.tocenter";
const string PhaseToCenterProcessor::NAME = "xform.phaseorigin.tocenter";
const string PhaseToCornerProcessor::NAME = "xform.phaseorigin.tocorner";
const string AutoMask2DProcessor::NAME = "mask.auto2d";
const string AutoMaskAsymUnit::NAME = "mask.asymunit";
const string AutoMask3DProcessor::NAME = "mask.auto3d.thresh";
const string AutoMask3D2Processor::NAME = "mask.auto3d";
const string AutoMaskDustProcessor::NAME = "mask.dust3d";
const string AddMaskShellProcessor::NAME = "mask.addshells";
const string IterMultiMaskProcessor::NAME = "mask.addshells.multilevel";
const string IterBinMaskProcessor::NAME = "mask.addshells.gauss";
const string PhaseToMassCenterProcessor::NAME = "xform.phasecenterofmass";
const string ToMassCenterProcessor::NAME = "xform.centerofmass";
const string ToCenterProcessor::NAME = "xform.center";
const string ACFCenterProcessor::NAME = "xform.centeracf";
const string SNRProcessor::NAME = "eman1.filter.snr";
const string FileFourierProcessor::NAME = "eman1.filter.byfile";
const string FSCFourierProcessor::NAME = "filter.wiener.byfsc";
const string SymSearchProcessor::NAME = "misc.symsearch";
const string LocalNormProcessor::NAME = "normalize.local";
const string StripeXYProcessor::NAME = "math.xystripefix";
const string BadLineXYProcessor::NAME = "math.xybadline";
const string IndexMaskFileProcessor::NAME = "mask.fromfile";
const string CoordinateMaskFileProcessor::NAME = "mask.fromfile.sizediff";
const string PaintProcessor::NAME = "mask.paint";
const string DirectionalSumProcessor::NAME = "misc.directional_sum";
template<> const string BinaryOperateProcessor<MaxPixelOperator>::NAME = "math.max";		// These 2 should not really be processors
template<> const string BinaryOperateProcessor<MinPixelOperator>::NAME = "math.min";
const string MaxPixelOperator::NAME = "math.max";
const string MinPixelOperator::NAME = "math.min";
const string MatchSFProcessor::NAME = "filter.matchto";
const string SetSFProcessor::NAME = "filter.setstrucfac";
const string SmartMaskProcessor::NAME = "mask.smart";
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
const string MirrorProcessor::NAME = "xform.mirror";
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
const string LowpassRandomPhaseProcessor::NAME = "filter.lowpass.randomphase";
const string NewLowpassButterworthProcessor::NAME = "filter.lowpass.butterworth";
const string NewHighpassButterworthProcessor::NAME = "filter.highpass.butterworth";
const string NewHomomorphicButterworthProcessor::NAME = "filter.homomorphic.butterworth";
const string NewLowpassTanhProcessor::NAME = "filter.lowpass.tanh";
const string NewHighpassTanhProcessor::NAME = "filter.highpass.tanh";
const string NewHomomorphicTanhProcessor::NAME = "filter.homomorphic.tanh";
const string NewBandpassTanhProcessor::NAME = "filter.bandpass.tanh";
const string CTF_Processor::NAME = "filter.CTF_";
const string ConvolutionKernelProcessor::NAME = "filter.convolution.kernel";
const string RotateInFSProcessor::NAME = "rotateinfs";
const string CircularAverageBinarizeProcessor::NAME = "threshold.binary.circularmean";
const string ObjDensityProcessor::NAME = "morph.object.density";
const string ObjLabelProcessor::NAME = "morph.object.label";
const string BwThinningProcessor::NAME = "morph.thin";
const string BwMajorityProcessor::NAME = "morph.majority";
const string PruneSkeletonProcessor::NAME = "morph.prune"; 
//#ifdef EMAN2_USING_CUDA
//const string CudaMultProcessor::NAME = "cuda.math.mult";
//const string CudaCorrelationProcessor::NAME = "cuda.correlate";
//#endif //EMAN2_USING_CUDA

#if 0
//const string XYZProcessor::NAME = "XYZ";
#endif	//0


template <> Factory < Processor >::Factory()
{
	force_add<HighpassAutoPeakProcessor>();
	force_add<LinearRampFourierProcessor>();
	force_add<LoGFourierProcessor>();
	force_add<DoGFourierProcessor>();
	force_add<AzSharpProcessor>();

	force_add<AmpweightFourierProcessor>();
	force_add<Axis0FourierProcessor>();
	force_add<Wiener2DFourierProcessor>();
	force_add<LowpassAutoBProcessor>();

	force_add<LinearPyramidProcessor>();
	force_add<LinearRampProcessor>();
	force_add<AbsoluateValueProcessor>();
	force_add<FloorValueProcessor>();
	force_add<BooleanProcessor>();
	force_add<KmeansSegmentProcessor>();
	force_add<DistanceSegmentProcessor>();
	force_add<ValuePowProcessor>();
	force_add<ValueSquaredProcessor>();
	force_add<ValueSqrtProcessor>();
	force_add<Rotate180Processor>();
	force_add<TransformProcessor>();
	force_add<ScaleTransformProcessor>();
	force_add<ApplySymProcessor>();
	force_add<IntTranslateProcessor>();
	force_add<RecipCarefullyProcessor>();
	force_add<SubtractOptProcessor>();

	force_add<ClampingProcessor>();
	force_add<NSigmaClampingProcessor>();

	force_add<ToZeroProcessor>();
	force_add<AboveToZeroProcessor>();
	force_add<OutlierProcessor>();
	force_add<ToMinvalProcessor>();
	force_add<CutToZeroProcessor>();
	force_add<BinarizeProcessor>();
//	force_add<BinarizeAmpProcessor>();
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
	force_add<MaskSoftProcessor>();
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
	force_add<NonConvexProcessor>();

	force_add<MakeRadiusSquaredProcessor>();
	force_add<MakeRadiusProcessor>();

	force_add<ComplexNormPixel>();

	force_add<LaplacianProcessor>();
//	force_add<ZeroConstantProcessor>();   // this is badly written, it does not work and never did. Who wrote this !?!?

	force_add<BoxMedianProcessor>();
	force_add<BoxSigmaProcessor>();
	force_add<BoxMaxProcessor>();

	force_add<MinusPeakProcessor>();
	force_add<PeakOnlyProcessor>();
	force_add<DiffBlockProcessor>();

	force_add<CutoffBlockProcessor>();
//	force_add<GradientRemoverProcessor>();
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
	force_add<AutoMaskDustProcessor>();
	force_add<AddMaskShellProcessor>();
	force_add<IterMultiMaskProcessor>();
	force_add<IterBinMaskProcessor>();
	force_add<AutoMaskAsymUnit>();

	force_add<CTFSNRWeightProcessor>();

	force_add<ToMassCenterProcessor>();
	force_add<ToCenterProcessor>();
	force_add<PhaseToMassCenterProcessor>();
	force_add<ACFCenterProcessor>();
//	force_add<SNRProcessor>();
	force_add<FSCFourierProcessor>();

	force_add<XGradientProcessor>();
	force_add<YGradientProcessor>();
	force_add<ZGradientProcessor>();

//	force_add<FileFourierProcessor>();

	force_add<SymSearchProcessor>();
	force_add<StripeXYProcessor>();
	force_add<BadLineXYProcessor>();
	force_add<LocalNormProcessor>();

	force_add<IndexMaskFileProcessor>();
	force_add<CoordinateMaskFileProcessor>();
	force_add<SetSFProcessor>();
	force_add<MatchSFProcessor>();

	force_add<SmartMaskProcessor>();

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
	force_add<LowpassRandomPhaseProcessor>();
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
	force_add<ConvolutionKernelProcessor>();

	//Gorgon-related processors
	force_add<ModelEMCylinderProcessor>();
	force_add<ApplyPolynomialProfileToHelix>();
	force_add<BinarySkeletonizerProcessor>();
	force_add<RotateInFSProcessor>();
	force_add<CircularAverageBinarizeProcessor>();
	force_add<ObjDensityProcessor>();
	force_add<ObjLabelProcessor>();
	force_add<BwThinningProcessor>();
	force_add<BwMajorityProcessor>();
	force_add<PruneSkeletonProcessor>();


//#ifdef EMAN2_USING_CUDA
//	force_add<CudaMultProcessor>();
//	force_add<CudaCorrelationProcessor>();
//#endif // EMAN2_USING_CUDA

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
//	printf("Default copy called\n");
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

	bool return_radial=(bool)params.set_default("return_radial",0);
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
	if (return_radial) image->set_attr("filter_curve",yarray);

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

	bool return_radial=(bool)params.set_default("return_radial",0);
	bool interpolate=(bool)params.set_default("interpolate",1);

	float cornerscale;
	if (image->get_zsize()>1) cornerscale=sqrt(3.0);
	else cornerscale=sqrt(2.0);
	
	if (image->is_complex()) {
		vector <float>yarray = image->calc_radial_dist(floor(image->get_ysize()*cornerscale/2),0,1.0,1);
		create_radial_func(yarray,image);
		image->apply_radial_func(0, 0.5f/image->get_ysize(), yarray,interpolate);
		if (return_radial) image->set_attr("filter_curve",yarray);
	}
	else {
		EMData *fft = image->do_fft();
		vector <float>yarray = fft->calc_radial_dist((int)floor(fft->get_ysize()*cornerscale/2.0),0,1.0,1);
		create_radial_func(yarray,fft);
		// 4/30/10 stevel turned off interpolation to fix problem with matched filter
		// 9/12/14 stevel, not sure why I turned off interp. Seems to cause rather than fix problems. Adding option to enable. Afraid to turn it on
		fft->apply_radial_func(0, 1.0f/image->get_ysize(), yarray,interpolate);
		EMData *ift = fft->do_ift();

		memcpy(image->get_data(),ift->get_data(),ift->get_xsize()*ift->get_ysize()*ift->get_zsize()*sizeof(float));
		if (return_radial) image->set_attr("filter_curve",yarray);

//		for (int i=0; i<yarray.size(); i++) printf("%d\t%f\n",i,yarray[i]);
		
		//ift->update(); Unecessary

		delete fft;
		delete ift;

	}

	image->update();
}

void AzSharpProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	
	float az_scale=(float)params.set_default("az_scale",1.0);
	float cornerscale;
	EMData *fft;
	
	if (image->is_complex()) fft=image;
	else EMData *fft = image->do_fft();

	int nx=fft->get_xsize();
	int ny=fft->get_ysize();
	int nz=fft->get_zsize();
	if (nz!=1) {
	}
	else {
		
	}
	if (image->is_complex()) {
		EMData *ift = fft->do_ift();

		memcpy(image->get_data(),ift->get_data(),ift->get_xsize()*ift->get_ysize()*ift->get_zsize()*sizeof(float));

		delete fft;
		delete ift;
	}

	image->update();
}


void LowpassAutoBProcessor::create_radial_func(vector < float >&radial_mask,EMData *image) const{
	float apix=(float)image->get_attr("apix_x");
	int verbose=(int)params["verbose"];
//	int adaptnoise=params.set_default("adaptnoise",0);
	float noisecutoff=(float)params.set_default("noisecutoff",1.0/6.0);
	if (apix<=0 || apix>7.0f) throw ImageFormatException("apix_x > 7.0 or <0");
	float ds=1.0f/(apix*image->get_ysize());	// 0.5 is because radial mask is 2x oversampled
	unsigned int start=(int)floor(1.0/(15.0*ds));
	unsigned int end=radial_mask.size()-2;
	if (noisecutoff>0) end=(int)floor(noisecutoff/ds);
	if (end>radial_mask.size()-2) {
		if (verbose) printf("WARNING: specified noisecutoff too close to Nyquist, reset !");
		end=radial_mask.size()-2;
	}
	if (end<start+2) {
		printf("WARNING: noise cutoff too close to 15 A ! Results will not be good...");
		start=end-5;
	}

	FILE *out=NULL;
	if (verbose>2)  {
		printf("Autob -> %d - %d  ds=%g apix=%g rdlmsk=%d\n",start,end,ds,apix,int(radial_mask.size()));
		out=fopen("fitplot.txt","w");
	}
	int N=(radial_mask.size()-start-2);
	float *x=(float *)malloc(N*sizeof(float));
	float *y=(float *)malloc(N*sizeof(float));
	float *dy=(float *)malloc(N*sizeof(float));
	for (unsigned int i=start; i<radial_mask.size()-2; i++ ) {		// -2 is arbitrary because sometimes the last pixel or two have funny values
		x[i-start]=ds*ds*i*i;
		if (radial_mask[i]>0) y[i-start]=log(radial_mask[i]); // ok
		else if (i>start) y[i-start]=y[i-start-1];		// not good
		else y[i-start]=0.0;							// bad
		if (i>start &&i<radial_mask.size()-3) dy[i-start]=y[i-start]-y[i-start-1];	// creates a 'derivative' of sorts, for use in adaptnoise
		if (out) fprintf(out,"%f\t%f\n",x[i-start],y[i-start]);
	}
	if (out) fclose(out);

	float slope=0,intercept=0;
	Util::calc_least_square_fit(end-start, x,y,&slope,&intercept,1);

	if (verbose) printf("slope=%f  intercept=%f\n",slope,intercept);

	float B=(float)params["bfactor"]+slope*4.0f/2.0f;		// *4 is for Henderson definition, 2.0 is for intensity vs amplitude
	float B2=(float)params["bfactor"];

	if (verbose) printf("User B = %1.2f  Corrective B = %1.2f  Total B=%1.3f\n",(float)params["bfactor"],slope*2.0,B);

	float cutval=exp(-B*pow(end*ds,2.0f)/4.0f)/exp(-B2*pow(end*ds,2.0f)/4.0f);
	for (unsigned int i=0; i<radial_mask.size(); i++) {
		if (i<=end) radial_mask[i]=exp(-B*pow(i*ds,2.0f)/4.0f);
		else radial_mask[i]=cutval*exp(-B2*pow(i*ds,2.0f)/4.0f);
	}
	if (verbose>1) Util::save_data(0,ds,radial_mask,"filter.txt");

	free(x);
	free(y);
	free(dy);
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

void Axis0FourierProcessor::process_inplace(EMData * image)
{
	EMData *fft;
	float *fftd;
	int f=0;
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

	int nx=fft->get_xsize();
	int ny=fft->get_ysize();
	if (params.set_default("x",1)) {
		for (int x=2; x<nx; x++) fftd[x]=0;
	}
	if (params.set_default("y",1)) {
		for (int y=1; y<ny; y++) { fftd[y*nx]=0; fftd[y*nx+1]=0; }
	}

	if (f) {
		fft->update();
		EMData *ift=fft->do_ift();
		memcpy(image->get_data(),ift->get_data(),(nx-2)*ny*sizeof(float));
		delete fft;
		delete ift;
	}

	image->update();

}

void AmpweightFourierProcessor::process_inplace(EMData * image)
{
	EMData *fft;
	float *fftd;
	int f=0;
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
	size_t n = (size_t)fft->get_xsize()*fft->get_ysize()*fft->get_zsize();
	for (size_t i=0; i<n; i+=2) {
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

void DistanceSegmentProcessor::process_inplace(EMData *)
{
	printf("Process inplace not implemented. Please use process.\n");
	return;
}


EMData *DistanceSegmentProcessor::process(const EMData * const image)
{
	EMData * result = image->copy();

	float thr = params.set_default("thr",0.9f);
	float minsegsep = params.set_default("minsegsep",5.0f);
	float maxsegsep = params.set_default("maxsegsep",5.1f);
	int verbose = params.set_default("verbose",0);

	vector<Pixel> pixels=image->calc_highest_locations(thr);

	vector<float> centers(3);	// only 1 to start
	int nx=image->get_xsize();
	int ny=image->get_ysize();
	int nz=image->get_zsize();
//	int nxy=nx*ny;

	// seed the process with the highest valued point
	centers[0]=(float)pixels[0].x;
	centers[1]=(float)pixels[0].y;
	centers[2]=(float)pixels[0].z;
	pixels.erase(pixels.begin());

	// outer loop. We add one center per iteration
	// This is NOT a very efficient algorithm, it assumes points are fairly sparse
	while (pixels.size()>0) {
		// iterate over pixels until we find a new center (then stop), delete any 'bad' pixels
		// no iterators because we remove elements

		for (unsigned int i=0; i<pixels.size(); i++) {

			Pixel p=pixels[i];
			// iterate over existing centers to see if this pixel should be removed ... technically we only should need to check the last center
			for (unsigned int j=0; j<centers.size(); j+=3) {
				float d=Util::hypot3(centers[j]-p.x,centers[j+1]-p.y,centers[j+2]-p.z);
				if (d<minsegsep) {		// conflicts with existing center, erase
					pixels.erase(pixels.begin()+i);
					i--;
					break;
				}
			}
		}

		int found=0;
		for (unsigned int i=0; i<pixels.size() && found==0; i++) {
			Pixel p=pixels[i];

			// iterate again to see if this may be a new valid center. Start at the end so we tend to build chains
			for (unsigned int j=centers.size()-3; j>0; j-=3) {
				float d=Util::hypot3(centers[j]-p.x,centers[j+1]-p.y,centers[j+2]-p.z);
				if (d<maxsegsep) {		// we passed minsegsep question already, so we know we're in the 'good' range
					centers.push_back((float)p.x);
					centers.push_back((float)p.y);
					centers.push_back((float)p.z);
					pixels.erase(pixels.begin()+i);	// in the centers list now, don't need it any more
					found=1;
					break;
				}
			}
		}

		// If we went through the whole list and didn't find one, we need to reseed again
		if (!found && pixels.size()) {
			if (verbose) printf("New chain\n");
			centers.push_back((float)pixels[0].x);
			centers.push_back((float)pixels[0].y);
			centers.push_back((float)pixels[0].z);
			pixels.erase(pixels.begin());
		}

		if (verbose) printf("%d points found\n",(int)(centers.size()/3));
	}

	// after we have our list of centers classify pixels
	for (int z=0; z<nz; z++) {
		for (int y=0; y<ny; y++) {
			for (int x=0; x<nz; x++) {
				if (image->get_value_at(x,y,z)<thr) {
					result->set_value_at(x,y,z,-1.0);		//below threshold -> -1 (unclassified)
					continue;
				}
				int bcls=-1;			// best matching class
				float bdist=(float)(nx+ny+nz);	// distance for best class
				for (unsigned int c=0; c<centers.size()/3; c++) {
					float d=Util::hypot3(x-centers[c*3],y-centers[c*3+1],z-centers[c*3+2]);
					if (d<bdist) { bdist=d; bcls=c; }
				}
				result->set_value_at(x,y,z,(float)bcls);		// set the pixel to the class number
			}
		}
	}

	result->set_attr("segment_centers",centers);

	return result;
}

EMData* ApplySymProcessor::process(const EMData * const image)
{
	Averager* imgavg = Factory<Averager>::get((string)params.set_default("avger","mean"));

	if (image->get_zsize()==1) {
		string s=(string)params["sym"];
		if (s[0]!='c' && s[0]!='C') throw ImageDimensionException("xform.applysym: Cn symmetry required for 2-D symmetrization");
		int n=atoi(s.c_str()+1);
		if (n<=0) throw InvalidValueException(n,"xform.applysym: Cn symmetry, n>0");

		for (int i=0; i<n; i++) {
			Transform t(Dict("type","2d","alpha",(float)(i*360.0f/n)));
			EMData* transformed = image->process("xform",Dict("transform",&t));
			imgavg->add_image(transformed);
			delete transformed;
		}
		EMData *ret=imgavg->finish();
		delete imgavg;
		return ret;
	}

	Symmetry3D* sym = Factory<Symmetry3D>::get((string)params.set_default("sym","c1"));
	vector<Transform> transforms = sym->get_syms();

	for(vector<Transform>::const_iterator trans_it = transforms.begin(); trans_it != transforms.end(); trans_it++) {
		Transform t = *trans_it;
		EMData* transformed = image->process("xform",Dict("transform",&t));
		imgavg->add_image(transformed);
		delete transformed;
	}
	EMData *ret=imgavg->finish();
	delete imgavg;
	return ret;
}

void ApplySymProcessor::process_inplace(EMData* image)
{
	EMData *tmp=process(image);
	memcpy(image->get_data(),tmp->get_data(),(size_t)image->get_xsize()*image->get_ysize()*image->get_zsize()*sizeof(float));
	delete tmp;
	image->update();
	return;
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
	bool pseudoatom = params.set_default("pseudoatom",0);
	float sep = params.set_default("sep",3.78f);

	int nx=image->get_xsize();
	int ny=image->get_ysize();
	int nz=image->get_zsize();
//	int nxy=nx*ny;

	// seed
	vector<float> centers(nseg*3);
	vector<float> count(nseg);
	// Alternative seeding method for paudoatom generation. Seed on the gird.
	if (pseudoatom){
		float ax=image->get_attr("apix_x");
		sep/=ax;
		if (verbose) printf("Seeding .....\n");
		int sx=int(nx/sep)+1,sy=int(ny/sep)+1,sz=int(nz/sep)+1;
		for(int setthr=0; setthr<10; setthr++){
			EMData m(sx,sy,sz);
			EMData mcount(sx,sy,sz);
			for (int i=0; i<nx; i++){
				for (int j=0; j<ny; j++){
					for (int k=0; k<nz; k++){
						int ni=(i/sep),nj=(j/sep),nk=(k/sep);
						float v=image->get_value_at(i,j,k);
						if (v>thr){
							m.set_value_at(ni,nj,nk,(m.get_value_at(ni,nj,nk)+v));
							mcount.set_value_at(ni,nj,nk,(mcount.get_value_at(ni,nj,nk)+1));
						}
					}
				}
			}
			m.div((nx/sx)*(ny/sy)*(nz/sz));
			int nsum=0;
			float l=image->get_attr("minimum"),r=image->get_attr("maximum"),th=0;
			while (abs(nsum-nseg)>0){
				th=(l+r)/2;
				nsum=0;
				for (int i=0; i<sx; i++){
					for (int j=0; j<sy; j++){
						for (int k=0; k<sz; k++){
							if (m.get_value_at(i,j,k)>th)  nsum+=1;
						}
					}
				}
				if (verbose) printf("%3f\t %3f\t %3f,\t %4d\t %4d\n", l,th,r,nsum,nseg);
				if (nsum>nseg) l=th;
				if (nsum<nseg) r=th;
				if ((r-l)<.01) break;
			}
	// 		nseg=nsum;
			int q=0;
			for (int i=0; i<sx; i++){
				for (int j=0; j<sy; j++){
					for (int k=0; k<sz; k++){
						if (m.get_value_at(i,j,k)>th){
							if(q<nseg*3){
								centers[q]=  float(i+.5)*sep;
								centers[q+1]=float(j+.5)*sep;
								centers[q+2]=float(k+.5)*sep;
								q+=3;
							}
						}
					}
				}
			}
			if (thr==-1.0e30f){
				printf("Estimated map threshold: %4f\n", th);
				thr=th;
			}
			else
				break;
		}
	}
	// Default: random seeding.
	else{
		for (int i=0; i<nseg*3; i+=3) {
			centers[i]=  Util::get_frand(0.0f,(float)nx);
			centers[i+1]=Util::get_frand(0.0f,(float)ny);
			centers[i+2]=Util::get_frand(0.0f,(float)nz);
		}
	}

	for (int iter=0; iter<maxiter; iter++) {
		// **** classify
		size_t pixmov=0;		// count of moved pixels
		for (int z=0; z<nz; z++) {
			for (int y=0; y<ny; y++) {
				for (int x=0; x<nx; x++) {
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
				for (int x=0; x<nx; x++) {
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

void KmeansSegmentProcessor::process_inplace(EMData *)
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
// 	float *fftd;
// 	int f=0;

	if (!image) {
		LOGWARN("NULL Image");
		return ret;
	}
	throw NullPointerException("Processor not yet implemented");

// 	if (!image->is_complex()) {
// 		fft = image->do_fft();
// 		fftd = fft->get_data();
// 		f=1;
// 	}
// 	else {
// 		fft=image;
// 		fftd=image->get_data();
// 	}

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

void LowpassRandomPhaseProcessor::create_radial_func(vector < float >&radial_mask) const { };

void LowpassRandomPhaseProcessor::process_inplace(EMData *image)
{
	float cutoff=0;
	preprocess(image);
	if( params.has_key("cutoff_abs") ) {
		cutoff=(float)params["cutoff_abs"];
	}
	else {
		printf("A cutoff_* parameter is required by filter.lowpass.randomphase\n");
		return;
	}


	if (image->get_zsize()==1) {
		int flag=0;
		if (!image->is_complex()) { image->do_fft_inplace(); flag=1; }
		image->ri2ap();
		int nx=image->get_xsize();
		int ny=image->get_ysize();

		int z=0;
		float *data=image->get_data();
		for (int y=-ny/2; y<ny/2; y++) {
			for (int x=0; x<nx/2+1; x++) {
				if (hypot(x/float(nx),y/float(ny))>=cutoff) {
					size_t idx=image->get_complex_index_fast(x,y,z);		// location of the amplitude
					data[idx+1]=Util::get_frand(0.0f,(float)(M_PI*2.0));
				}
			}
		}

		image->ap2ri();

		if (flag) {
			image->do_ift_inplace();
			image->depad();
		}
	}
	else {		// 3D
		int flag=0;
		if (!image->is_complex()) { image->do_fft_inplace(); flag=1; }
		image->ri2ap();
		int nx=image->get_xsize();
		int ny=image->get_ysize();
		int nz=image->get_zsize();

		float *data=image->get_data();
		for (int z=-nz/2; z<nz/2; z++) {
			for (int y=-ny/2; y<ny/2; y++) {
				for (int x=0; x<nx/2; x++) {
					if (Util::hypot3(x/float(nx),y/float(ny),z/float(nz))>=cutoff) {
						size_t idx=image->get_complex_index_fast(x,y,z);		// location of the amplitude
						data[idx+1]=Util::get_frand(0.0f,(float)(M_PI*2.0));
					}
				}
			}
		}
		image->ap2ri();

		if (flag) {
			image->do_ift_inplace();
			image->depad();
		}
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
		highpass = (float)params["cutoff_freq"] * (float)dict["apix_x"] * (float)dict["ny"] / 2.0f;
	}
	else if( params.has_key("cutoff_pixels") ) {
		highpass = (float)params["cutoff_pixels"] / (float)dict["nx"];
	}
}

void HighpassAutoPeakProcessor::create_radial_func(vector < float >&radial_mask, EMData *) const
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

void LoGFourierProcessor::create_radial_func(vector < float >&radial_mask) const
{

	Assert(radial_mask.size() > 0);
	float x = 0.0f , nqstep = 0.5f/radial_mask.size();
	size_t size=radial_mask.size();
	float var = sigma*sigma;
	for (size_t i = 0; i < size; i++) {
		radial_mask[i] = ((x*x - var)/var*var)*exp(-x*x/2*var);
		x += nqstep;
	}
}

void DoGFourierProcessor::create_radial_func(vector < float >&radial_mask) const
{

	Assert(radial_mask.size() > 0);
	float x = 0.0f , nqstep = 0.5f/radial_mask.size();
	size_t size=radial_mask.size();
	float norm = 1.0f/sqrt(2*M_PI);
	for (size_t i = 0; i < size; i++) {
		radial_mask[i] = norm*((1.0f/sigma1*exp(-x*x/(2.0f*sigma1*sigma1))) - (1.0f/sigma2*exp(-x*x/(2.0f*sigma2*sigma2))));
		x += nqstep;
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

	for (size_t i = 0; i < size; ++i) {
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



void CircularMaskProcessor::calc_locals(EMData *)
{
	xc = Util::fast_floor(nx/2.0f) + dx;
	yc = Util::fast_floor(ny/2.0f) + dy;
	zc = Util::fast_floor(nz/2.0f) + dz;

	if (outer_radius < 0) {
		outer_radius = nx / 2 + outer_radius +1;
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

	xc = Util::fast_floor(nx/2.0f) + dx;
	yc = Util::fast_floor(ny/2.0f) + dy;
	zc = Util::fast_floor(nz/2.0f) + dz;

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



	for (size_t i = 0; i < size; ++i) {
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

	size_t cpysize = areasize * sizeof(float);
	size_t start = (nx * ny + nx + 1) * n;

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

				k = (size_t)z * nsec + y * nx + x;

				for (int bz = zbox_start; bz < zbox_end; bz++) {
					for (int by = 0; by < areasize; by++) {
						memcpy(&matrix[(size_t)bz * box_nsec + by * areasize],
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
		memset(kernel, 0, (size_t)areasize * areasize * areasize);
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
		size_t knxy = (size_t)k * nxy;

		for (int j = n; j < ny - n; j++) {
			int jnx = j * nx;

			for (int i = n; i < nx - n; i++) {
				size_t s = 0;

				for (int i2 = i - n; i2 <= i + n; i2++) {
					for (int j2 = j - n; j2 <= j + n; j2++) {
						for (int k2 = k - nzz; k2 <= k + nzz; k2++) {
							array[s] = data2[i2 + j2 * nx + (size_t)k2 * nxy];
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
		size_t cur_l = (size_t)l * nxy_new;

		for (int j = 0; j < ny; j++) {
			int j_min = j * shrink_factor;
			int j_max = (j + 1) * shrink_factor;
			size_t cur_j = j * nx + cur_l;

			for (int i = 0; i < nx; i++) {
				int i_min = i * shrink_factor;
				int i_max = (i + 1) * shrink_factor;

				size_t k = 0;
				for (int l2 = l_min; l2 < l_max; l2++) {
					size_t cur_l2 = l2 * nxy_old;

					for (int j2 = j_min; j2 < j_max; j2++) {
						size_t cur_j2 = j2 * nx_old + cur_l2;

						for (int i2 = i_min; i2 < i_max; i2++) {
							mbuf[k] = data_copy[i2 + cur_j2];
							++k;
						}
					}
				}

				for (k = 0; k < size_t(threed_shrink_factor / 2 + 1); k++) {
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

	size_t nx = from->get_xsize();
	size_t ny = from->get_ysize();
	size_t nz = from->get_zsize();
	size_t nxy = nx*ny;


	size_t shrunken_nx = nx / shrink_factor;
	size_t shrunken_ny = ny / shrink_factor;
	size_t shrunken_nz = 1;
	size_t shrunken_nxy = shrunken_nx * shrunken_ny;

	int normalize_shrink_factor = shrink_factor * shrink_factor;
	int z_shrink_factor = 1;

	if (nz > 1) {
		shrunken_nz = nz / shrink_factor;
		normalize_shrink_factor *= shrink_factor;
		z_shrink_factor = shrink_factor;
	}

	float invnormfactor = 1.0f/(float)normalize_shrink_factor;

	for (size_t k = 0; k < shrunken_nz; k++) {
		size_t k_min = k * shrink_factor;
		size_t k_max = k * shrink_factor + z_shrink_factor;
		size_t cur_k = k * shrunken_nxy;

		for (size_t j = 0; j < shrunken_ny; j++) {
			size_t j_min = j * shrink_factor;
			size_t j_max = j * shrink_factor + shrink_factor;
			size_t cur_j = j * shrunken_nx + cur_k;

			for (size_t i = 0; i < shrunken_nx; i++) {
				size_t i_min = i * shrink_factor;
				size_t i_max = i * shrink_factor + shrink_factor;

				float sum = 0;
				for (size_t kk = k_min; kk < k_max; kk++) {
					size_t cur_kk = kk * nxy;

					for (size_t jj = j_min; jj < j_max; jj++) {
						size_t cur_jj = jj * nx + cur_kk;
						for (size_t ii = i_min; ii < i_max; ii++) {
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

void NonConvexProcessor::process_inplace(EMData * image) {
	if (!image) { LOGWARN("NULL IMAGE"); return; }
	//int isinten=image->get_attr_default("is_intensity",0);

	// 1-D
	if (image->get_ysize()==1) {

	}
	// 2-D
	else if (image->get_zsize()==1) {
//		if (!isinten) throw ImageDimensionException("Only complex intensity images currently supported by NonConvexProcessor");
		int nx2=image->get_xsize()/2;
		int ny2=image->get_ysize()/2;
		vector<float> rdist = image->calc_radial_dist(nx2*1.5,0,1,false);		// radial distribution to make sure nonconvex values decrease radially
		// Make sure rdist is decreasing (or flat)
		for (int i=1; i<nx2; i++) {
			if (rdist[i]>rdist[i-1]) rdist[i]=rdist[i-1];
		}

		image->process_inplace("xform.fourierorigin.tocenter");
		EMData* binary=image->copy();

		// First we eliminate convex points from the input image (set to zero)
		for (int x=0; x<image->get_xsize(); x+=2) {
			for (int y=1; y<image->get_ysize()-1; y++) {
				int r=(int)hypot((float)(x/2),(float)(y-ny2));
				float cen=(*binary)(x,y);
				if (x==0 || x==nx2*2-2 || (cen>(*binary)(x+2,y) || cen>(*binary)(x-2,y) || cen>(*binary)(x,y+1) || cen >(*binary)(x,y-1) || (*binary)(x,y)>rdist[r])) {		// point is considered nonconvex if lower than surrounding values and lower than mean
					image->set_value_at_fast(x/2+nx2,y,0.0);	// we are turning image into a full real-space intensity image for now
					image->set_value_at_fast(nx2-x/2,ny2*2-y-1,0.0);
				}
				else {
					image->set_value_at_fast(x/2+nx2,y,cen);	// we are turning image into a full real-space intensity image for now
					image->set_value_at_fast(nx2-x/2,ny2*2-y-1,cen);	// It will contain non-zero values only for nonconvex points
				}
			}
		}
		image->set_value_at_fast(nx2+1,ny2,(*binary)(2,ny2));	// We keep the points near the Fourier origin as a central anchor even though it's convex
		image->set_value_at_fast(nx2-1,ny2,(*binary)(2,ny2));	// We keep the points near the Fourier origin as a central anchor even though it's convex
		image->set_value_at_fast(nx2,ny2+1,(*binary)(0,ny2+1));	// We keep the points near the Fourier origin as a central anchor even though it's convex
		image->set_value_at_fast(nx2,ny2-1,(*binary)(0,ny2-1));	// We keep the points near the Fourier origin as a central anchor even though it's convex
		for (int y=0; y<ny2*2; y++) image->set_value_at_fast(0,y,0.0f);

		// Now make a binary version of the convex points
		float *idat=image->get_data();
		float *bdat=binary->get_data();
		int nxy=(nx2*ny2*4);
		for (int i=0; i<nxy; i++) {
			bdat[i]=idat[i]==0?0:1.0f;		// binary version of the convex points in image
		}
		binary->update();

		// We now use a Gaussian filter on both images, to use Gaussian interpolation to fill in zero values
		image->set_complex(false);		// so we can use a Gaussian filter on it
		binary->set_complex(false);

/*		image->write_image("con.hdf",0);*/
		image->set_fftpad(false);
		binary->set_fftpad(false);

		// Gaussian blur of both images
		image->process_inplace("filter.lowpass.gauss",Dict("cutoff_abs",0.04f));
		binary->process_inplace("filter.lowpass.gauss",Dict("cutoff_abs",0.04f));

/*		image->write_image("con.hdf",1);
		binary->write_image("con.hdf",2);*/

		for (int x=0; x<image->get_xsize(); x+=2) {
			for (int y=0; y<image->get_ysize(); y++) {
				float bv=binary->get_value_at(x/2+nx2,y);
				image->set_value_at_fast(x,y,image->get_value_at(x/2+nx2,y)/(bv<=0?1.0f:bv));
				image->set_value_at_fast(x+1,y,0.0);
			}
		}
		image->set_complex(true);
		image->set_fftpad(true);
		image->process_inplace("xform.fourierorigin.tocorner");
		delete binary;
	}
	else throw ImageDimensionException("3D maps not yet supported by NonConvexProcessor");

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

	float nonzero = params.set_default("nonzero",false);

	float *d = image->get_data();
	int i = 0;
	int j = 0;

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	float zval=9.99e23f;		// we're assuming we won't run into this exact value for an edge, not great programming, but good enough
	if (nonzero) {
		int x,y;
		size_t corn=nx*ny-1;

		// this set of 4 tests looks for any edges with exactly the same value
		for (x=1; x<nx; x++) { if (d[x]!=d[0]) break;}
		if (x==nx) zval=d[0];

		for (y=1; y<ny; y++) { if (d[y*nx]!=d[0]) break; }
		if (y==ny) zval=d[0];

		for (x=1; x<nx; x++) { if (d[corn-x]!=d[corn]) break;}
		if (x==nx) zval=d[corn];

		for (y=1; y<ny; y++) { if (d[corn-y*nx]!=d[corn]) break; }
		if (y==ny) zval=d[corn];

		if (zval!=9.99e23f) { image->set_attr("hadzeroedge",1); printf("zeroedge %f\n",zval); }
		else image->set_attr("hadzeroedge",0);

		// This tries to detect images where the edges have been filled with the nearest non-zero value. The filter does nothing, but we set the tag.
		for (x=nx/2-5; x<nx/2+5; x++) {
			if (d[x]!=d[x+nx] || d[x]!=d[x+nx*2] ) break;
		}
		if (x==nx/2+5) image->set_attr("hadzeroedge",2);

		for (x=nx/2-5; x<nx/2+5; x++) {
			if (d[corn-x]!=d[corn-x-nx] || d[corn-x]!=d[corn-x-nx*2]) break;
		}
		if (x==nx/2+5) image->set_attr("hadzeroedge",2);

		for (y=ny/2-5; y<ny/2+5; y++) {
			if (d[y*nx]!=d[y*nx+1] || d[y*nx]!=d[y*nx+2] ) break;
		}
		if (y==ny/2+5) image->set_attr("hadzeroedge",2);

		for (y=ny/2-5; y<ny/2+5; y++) {
			if (d[corn-y*nx]!=d[corn-y*nx-1] || d[corn-y*nx]!=d[corn-y*nx-2]) break;
		}
		if (y==ny/2+5) image->set_attr("hadzeroedge",2);

	}
	if (zval==9.99e23f) zval=0;

	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx - 1; i++) {
			if (d[i + j * nx] != zval) {
				break;
			}
		}

		float v = d[i + j * nx];
		while (i >= 0) {
			d[i + j * nx] = v;
			i--;
		}

		for (i = nx - 1; i > 0; i--) {
			if (d[i + j * nx] != zval)
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
			if (d[i + j * nx] != zval)
				break;
		}

		float v = d[i + j * nx];
		while (j >= 0) {
			d[i + j * nx] = v;
			j--;
		}

		for (j = ny - 1; j > 0; j--) {
			if (d[i + j * nx] != zval)
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

void OutlierProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	bool fix_zero=(bool)params.set_default("fix_zero",0);
	float sigmamult=(float)params.set_default("sigma",3.0);
	
	if (sigmamult<=0.0) throw InvalidValueException(sigmamult,"threshold.outlier.localmean: sigma must be >0");
	
	float hithr=(float)image->get_attr("mean")+(float)(image->get_attr("sigma"))*sigmamult;
	float lothr=(float)image->get_attr("mean")-(float)(image->get_attr("sigma"))*sigmamult;

	int nx=image->get_xsize();
	int ny=image->get_ysize();
	int nz=image->get_zsize();
	
	// This isn't optimally efficient
	EMData *im[2];
	im[0]=image;
	im[1]=image->copy_head();
	
	if (nz==1) {
		int repeat=1;
		while (repeat) {
			memcpy(im[1]->get_data(),im[0]->get_data(),image->get_xsize()*image->get_ysize()*image->get_zsize()*sizeof(float));
			repeat=0;
			for (int y=0; y<ny; y++) {
				for (int x=0; x<nx; x++) {
					// if the pixel is an outlier
					float pix=im[1]->get_value_at(x,y);
					if (pix>hithr || pix<lothr || (pix==0 && fix_zero)) {
						int y0=0>y-1?0:y-1;
						int y1=y>=ny-1?ny-1:y+1;
						int x0=0>x-1?0:x-1;
						int x1=x>=nx-1?nx-1:x+1;
						float c=0.0f,nc=0.0f;
						for (int yy=y0; yy<=y1; yy++) {
							for (int xx=x0; xx<=x1; xx++) {
								float lpix=im[1]->get_value_at(xx,yy);
								if (lpix>hithr || lpix<lothr || (lpix==0 && fix_zero)) continue;
								c+=lpix;
								nc++;
							}
						}
						if (nc!=0) im[0]->set_value_at(x,y,c/nc);
						else repeat=1;
					}
				}
			}
		}
	}
	else {
		throw ImageDimensionException("threshold.outlier.localmean: 3D not yet implemented");
	}
	delete im[1];
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

	size_t size = (size_t)image->get_xsize() * image->get_ysize() * image->get_zsize();
	float *data = image->get_data();

	for (size_t i = 0; i < size; ++i) {
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
		size_t size = (size_t)image->get_xsize() * image->get_ysize() * image->get_zsize();
		double sum = 0;
		double sq2 = 0;
		size_t n_norm = 0;

		for (size_t i = 0; i < size; ++i) {
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
	size_t size = (size_t)image->get_xsize() * image->get_ysize() * image->get_zsize();
	double sum = 0;
	size_t n_norm = 0;

	for (size_t i = 0; i < size; ++i) {
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
	int verbose = params.set_default("verbose",0);

	if (mass <= 0) throw InvalidParameterException("You must specify a positive non zero mass");

	float tthr = params.set_default("thr",(float)image->get_attr("mean")+(float)image->get_attr("sigma"));

	float apix = image->get_attr_default("apix_x",1.0f);
	apix = params.set_default("apix",apix);

	if (apix <= 0) throw InvalidParameterException("You must specify a positive non zero apix");

	float step = ((float)image->get_attr("sigma"))/5.0f;

	if (step==0) throw InvalidParameterException("This image has sigma=0, cannot give it mass");


	size_t n = image->get_size();
	float* d = image->get_data();


	float thr=(float)image->get_attr("mean")+(float)image->get_attr("sigma")/2.0;
	int count=0;
	for (size_t i=0; i<n; ++i) {
		if (d[i]>=thr) ++count;
	}
	if (verbose) printf("apix=%1.3f\tmass=%1.1f\tthr=%1.2f\tstep=%1.3g\n",apix,mass,thr,step);

	float max = image->get_attr("maximum");
	float min = image->get_attr("minimum");
	for (int j=0; j<4; j++) {
		int err=0;
		while (thr<max && count*apix*apix*apix*.81/1000.0>mass) {
			thr+=step;
			count=0;
			for (size_t i=0; i<n; ++i) {
				if (d[i]>=thr) ++count;
			}
			err+=1;
			if (err>1000) throw InvalidParameterException("Specified mass could not be achieved");
			if (verbose>1) printf("%d\t%d\t%1.3f\t%1.2f\n",err,count,thr,count*apix*apix*apix*.81/1000.0);
		}

		step/=4.0;

		while (thr>min && count*apix*apix*apix*.81/1000.0<mass) {
			thr-=step;
			count=0;
			for (size_t i=0; i<n; ++i) {
				if (d[i]>=thr) ++count;
			}
			err+=1;
			if (err>1000) throw InvalidParameterException("Specified mass could not be achieved");
			if (verbose>1) printf("%d\t%d\t%1.3f\t%1.2f\n",err,count,thr,count*apix*apix*apix*.81/1000.0);

		}
		step/=4.0;
	}

	image->mult((float)tthr/thr);
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
//	return image->get_circle_mean();
	int nx=image->get_xsize();
	int ny=image->get_ysize();
	int nz=image->get_zsize();

	float radius = params.set_default("radius",((float)ny/2-2));
	if (radius<0) radius=ny/2+radius;

 	static bool busy = false;		// avoid problems with threads and different image sizes
	static EMData *mask = 0;
	static int oldradius=radius;
	
	if (!mask || !EMUtil::is_same_size(image, mask)||radius!=oldradius) {
		while (busy) ;
		busy=true;
		if (!mask) {
			mask = new EMData();
		}
		mask->set_size(nx, ny, nz);
		mask->to_one();

		mask->process_inplace("mask.sharp", Dict("inner_radius", radius - 1,
							 "outer_radius", radius + 1));

	}
	busy=true;
	double n = 0,s=0;
	float *d = mask->get_data();
	float * data = image->get_data();
	size_t size = (size_t)nx*ny*nz;
	for (size_t i = 0; i < size; ++i) {
		if (d[i]) { n+=1.0; s+=data[i]; }
	}

	float result = (float)(s/n);
//	printf("cmean=%f\n",result);
	busy=false;
	
	return result;
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
	size_t nyz = ny * nz;

	for (size_t i = 0; i < nyz; i++) {
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

EMData *SubtractOptProcessor::process(const EMData * const image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return NULL;
	}

	EMData *refr = params["ref"];
	EMData *actual = params.set_default("actual",(EMData*)NULL);
	EMData *ref;
	bool return_radial = params.set_default("return_radial",false);
	bool return_fft = params.set_default("return_fft",false);
	bool ctfweight = params.set_default("ctfweight",false);
	bool return_presigma = params.set_default("return_presigma",false);
	bool return_subim = params.set_default("return_subim",false);
	int si0=(int)floor(params.set_default("low_cutoff_frequency",0.0f)*image->get_ysize());
	int si1=(int)ceil(params.set_default("high_cutoff_frequency",0.7071f)*image->get_ysize());		// include the corners unless explicitly excluded
	
	// We will be modifying imf, so it needs to be a copy
	EMData *imf;
	if (image->is_complex()) imf=image->copy();
	else imf=image->do_fft();

	if (ctfweight) {
		EMData *ctfi=imf->copy_head();
		Ctf *ctf;
//		if (image->has_attr("ctf")) 
			ctf=(Ctf *)(image->get_attr("ctf"));
//		else ctf=(Ctf *)(ref->get_attr("ctf"));
		ctf->compute_2d_complex(ctfi,Ctf::CTF_INTEN);
		imf->mult(*ctfi);
		delete ctfi;
	}
	
	// Make sure ref is complex
	if (refr->is_complex()) ref=refr;
	else ref=refr->do_fft();

	EMData *actf;
	if (actual==NULL) actf=ref;
	else {
		if (ctfweight) throw InvalidCallException("math.sub.optimal: Sorry, cannot use ctfweight in combination with actual");
		if (actual->is_complex()) actf=actual;
		else actf=actual->do_fft();
	}
		
	int ny2=(int)(image->get_ysize()*sqrt(2.0)/2);
	vector <double>rad(ny2+1);
	vector <double>norm(ny2+1);

	// We are essentially computing an FSC here, but while the reference (the image
	// we plan to subtract) is normalized, the other image is not. This gives us a filter
	// to apply to the reference to optimally eliminate its contents from 'image'.
	for (int y=-ny2; y<ny2; y++) {
		for (int x=0; x<ny2; x++) {
			int r=int(Util::hypot_fast(x,y));
			if (r>ny2) continue;
			std::complex<float> v1=imf->get_complex_at(x,y);
			std::complex<float> v2=ref->get_complex_at(x,y);
			rad[r]+=(double)(v1.real()*v2.real()+v1.imag()*v2.imag());
//			norm[r]+=v2.real()*v2.real()+v2.imag()*v2.imag()+v1.real()*v1.real()+v1.imag()*v1.imag();
			norm[r]+=(double)(v2.real()*v2.real()+v2.imag()*v2.imag());
		}
	}
	for (int i=0; i<ny2; i++) rad[i]/=norm[i];

// 	FILE *out=fopen("dbug.txt","w");
// 	for (int i=0; i<ny2; i++) fprintf(out,"%lf\t%lf\t%lf\n",(float)i,rad[i],norm[i]);
// 	fclose(out);
	
	float oldsig=-1.0;
	// This option computes the real-space sigma on the input-image after the specified filter
	// This is an expensive option, but more efficient than computing the same using other means
	if (return_presigma) {
		for (int y=-ny2; y<ny2; y++) {
			for (int x=0; x<imf->get_xsize()/2; x++) {
				int r=int(Util::hypot_fast(x,y));
				if (r>=ny2 || r>=si1 || r<si0) {
					imf->set_complex_at(x,y,0);
					continue;
				}
				std::complex<float> v1=imf->get_complex_at(x,y);
				imf->set_complex_at(x,y,v1);
			}
		}
		EMData *tmp=imf->do_ift();
		oldsig=(float)tmp->get_attr("sigma");
		delete tmp;
	}
	
	for (int y=-ny2; y<ny2; y++) {
		for (int x=0; x<imf->get_xsize()/2; x++) {
			int r=int(Util::hypot_fast(x,y));
			if (r>=ny2 || r>=si1 || r<si0) {
				imf->set_complex_at(x,y,0);
				continue;
			}
			std::complex<float> v1=imf->get_complex_at(x,y);
			std::complex<float> v2=actf->get_complex_at(x,y);
			v2*=(float)rad[r];
			if (return_subim) imf->set_complex_at(x,y,v2);
			else imf->set_complex_at(x,y,v1-v2);
		}
	}
	
	if (!refr->is_complex()) delete ref;
	if (actual!=NULL && !actual->is_complex()) delete actf;

	vector <float>radf;
	if (return_radial) {
		radf.resize(ny2);
		for (int i=0; i<ny2; i++) radf[i]=(float)rad[i];
	}
		
	if (!return_fft) {
		EMData *ret=imf->do_ift();
		delete imf;
		if (return_radial) ret->set_attr("filter_curve",radf);
		if (return_presigma) {
			ret->set_attr("sigma_presub",oldsig);
//			printf("Set %f\n",(float)image->get_attr("sigma_presub"));
		}
		return ret;
	}
	if (return_radial) imf->set_attr("filter_curve",radf);
	if (return_presigma) imf->set_attr("sigma_presub",oldsig);
	return imf;
}

void SubtractOptProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL image");
		return;
	}

	EMData *tmp=process(image);
	memcpy(image->get_data(),tmp->get_data(),(size_t)image->get_xsize()*image->get_ysize()*image->get_zsize()*sizeof(float));
	delete tmp;
	image->update();
	return;
}


void NormalizeToLeastSquareProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	EMData *to = params["to"];

	bool ignore_zero = params.set_default("ignore_zero",true);
	float ignore_lowsig = params.set_default("ignore_lowsig",-1.0);
	float low_threshold = params.set_default("low_threshold",-FLT_MAX);
	float high_threshold = params.set_default("high_threshold",FLT_MAX);

	float *dimage = image->get_data();
	float *dto = to->get_data();

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	size_t size = (size_t)nx * ny * nz;

	// rewrote this to just use GSL and get rid of David's old code.
	// The two passes are just to make sure we don't eat too much RAM if we do 3D
	if (ignore_lowsig<0) ignore_lowsig=0;
//	FILE *dbug = fopen("dbug.txt","w");

	size_t count=0;
	float meani=(float)image->get_attr("mean");
	float meant=(float)to->get_attr("mean");
	float sigi=(float)image->get_attr("sigma")*ignore_lowsig;
	float sigt=(float)to->get_attr("sigma")*ignore_lowsig;
	for (size_t i = 0; i < size; ++i) {
		if (dto[i] >= low_threshold && dto[i] <= high_threshold
			&& (dto[i]>=meant+sigt || dto[i]<=meant-sigt)
			&& (dimage[i]>=meani+sigi || dimage[i]<=meani-sigi)
			&& (!ignore_zero ||(dto[i] != 0.0f && dimage[i] != 0.0f))) {
			count++;
//			fprintf(dbug,"%f\t%f\n",dimage[i],dto[i]);
		}
	}
//	fclose(dbug);

	double *x=(double *)malloc(count*sizeof(double));
	double *y=(double *)malloc(count*sizeof(double));
	count=0;
	for (size_t i = 0; i < size; ++i) {
		if (dto[i] >= low_threshold && dto[i] <= high_threshold
			&& (dto[i]>=meant+sigt || dto[i]<=meant-sigt)
			&& (dimage[i]>=meani+sigi || dimage[i]<=meani-sigi)
			&& (!ignore_zero ||(dto[i] != 0.0f && dimage[i] != 0.0f))) {
			x[count]=dimage[i];
			y[count]=dto[i];
			count++;
		}
	}
	double c0,c1;
	double cov00,cov01,cov11,sumsq;
	gsl_fit_linear (x, 1, y, 1, count, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

	free(x);
	free(y);

	for (size_t i = 0; i < size; ++i) dimage[i]=dimage[i]*c1+c0;
	image->set_attr("norm_mult",c1);
	image->set_attr("norm_add",c0);
	image->update();
}

void BinarizeFourierProcessor::process_inplace(EMData* image) {
	ENTERFUNC;
	if (!image->is_complex()) throw ImageFormatException("Fourier binary thresholding processor only works for complex images");

	float threshold = params.set_default("value",0.0f);
	image->ri2ap(); //  works for cuda

	float* d = image->get_data();
	for( size_t i = 0; i < image->get_size()/2; ++i, d+=2) {
		if ( *d < threshold ) {
			*d = 0;
			*(d+1) = 0;
		}
	}
        image->ap2ri();
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
				size_t cur_k1 = (size_t)(k + half_width) * new_slice_size * is_3d;
				int cur_k2 = k * slice_size;

				for (int i = 0; i < height; i++) {
					int cur_i1 = (i + half_width) * new_width;
					int cur_i2 = i * width;

					for (int j = 0; j < width; j++) {
						size_t k1 = cur_k1 + cur_i1 + (j + half_width);
						int k2 = cur_k2 + cur_i2 + j;
						old_img[k1] = new_img[k2];
					}
				}
			}

			for (int k = 0; k < slicenum; k++) {
				size_t cur_k = (k + half_width) * new_slice_size * is_3d;

				for (int i = 0; i < height; i++) {
					int cur_i = (i + half_width) * new_width;

					for (int j = 0; j < half_width; j++) {
						size_t k1 = cur_k + cur_i + j;
						size_t k2 = cur_k + cur_i + (2 * half_width - j);
						old_img[k1] = old_img[k2];
					}

					for (int j = 0; j < half_width; j++) {
						size_t k1 = cur_k + cur_i + (width + half_width + j);
						size_t k2 = cur_k + cur_i + (width + half_width - j - 2);
						old_img[k1] = old_img[k2];
					}
				}


				for (int i = 0; i < half_width; i++) {
					int i2 = i * new_width;
					int i3 = (2 * half_width - i) * new_width;
					for (int j = 0; j < (width + 2 * half_width); j++) {
						size_t k1 = cur_k + i2 + j;
						size_t k2 = cur_k + i3 + j;
						old_img[k1] = old_img[k2];
					}

					i2 = (height + half_width + i) * new_width;
					i3 = (height + half_width - 2 - i) * new_width;
					for (int j = 0; j < (width + 2 * half_width); j++) {
						size_t k1 = cur_k + i2 + j;
						size_t k2 = cur_k + i3 + j;
						old_img[k1] = old_img[k2];
					}
				}
			}

			size_t idx;
			for (int k = 0; k < slicenum; k++) {
				size_t cur_k = (k + half_width) * new_slice_size;

				for (int i = 0; i < height; i++) {
					int cur_i = (i + half_width) * new_width;

					for (int j = 0; j < width; j++) {
						float f1 = 0;
						float f2 = 0;
						size_t k1 = cur_k + cur_i + (j + half_width);

						for (int p = zstart; p <= zend; p++) {
							size_t cur_p1 = (p + half_width) * (2 * half_width + 1) * (2 * half_width + 1);
							size_t cur_p2 = (k + half_width + p) * new_slice_size;

							for (int m = -half_width; m <= half_width; m++) {
								size_t cur_m1 = (m + half_width) * (2 * half_width + 1);
								size_t cur_m2 = cur_p2 + cur_i + m * new_width + j + half_width;

								for (int n = -half_width; n <= half_width; n++) {
									size_t k = cur_p1 + cur_m1 + (n + half_width);
									size_t k2 = cur_m2 + n;
									float f3 = Util::square(old_img[k1] - old_img[k2]);

									f3 = mask[k] * (1.0f / (1 + f3 / value_sigma));
									f1 += f3;
									size_t l1 = cur_m2 + n;
									f2 += f3 * old_img[l1];
								}

								idx = (size_t)k * height * width + i * width + j;
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
	if (EMData::usecuda == 1 && image->getcudarwdata()) {
		//cout << "flip processor" << endl;
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
				std::fill(d+z*nxy,d+(size_t)z*nxy+nx,0); // So if we have change it to the mean it's easy to do so. (instead of using memset)
			}
			for(int y=offset; y<ny/2; ++y) {
				for(int x=0; x<nx; ++x) {
					std::swap(d[(size_t)z*nxy + y*nx +x], d[(size_t)z*nxy + (ny -y -1+offset)*nx +x]);
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
					idx1 = (size_t)z*nxy + y*nx + x;
					idx2 = (size_t)(nz-z-1+offset)*nxy + y*nx + x;
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
				idx = (size_t)s*nxy+ny/2*nx+c;
				prev[0] = rdata[idx];
				prev[1] = rdata[idx+1];
				for( int r = 0; r <= ny/2; ++r ) {
					idx = (size_t)s*nxy+r*nx+c;
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
				idx1 = (size_t)s*nxy+r*nx+c;
				idx2 = (size_t)s*nxy+(r+ny/2)*nx+c;
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
					idx = (size_t)nz/2*nxy+r*nx+c;
					prev[0] = rdata[idx];
					prev[1] = rdata[idx+1];
					for( int s = 0; s <= nz/2; ++s ) {
						idx = (size_t)s*nxy+r*nx+c;
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
					idx1 = (size_t)s*nxy+r*nx+c;
					idx2 = (size_t)(s+nz/2)*nxy+r*nx+c;
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
						size_t i=x+y*nx+(size_t)z*nxy;
						size_t i2=x+y2*nx+(size_t)z2*nxy;
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
				idx = (size_t)s*nxy+c;
				prev[0] = rdata[idx];
				prev[1] = rdata[idx+1];
				for( int r = ny/2; r >= 0; --r ) {
					idx = (size_t)s*nxy+r*nx+c;
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
				idx1 = (size_t)s*nxy+r*nx+c;
				idx2 = (size_t)s*nxy+(r+ny/2+yodd)*nx+c;
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
						idx = (size_t)s*nxy+r*nx+c;
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
					idx1 = (size_t)s*nxy+r*nx+c;
					idx2 = (size_t)(s+nz/2+zodd)*nxy+r*nx+c;
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
		size_t k2 = (size_t)k * nxy;

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
					idx1 = (size_t)s*nxy+(size_t)r*nx+c;
					idx2 = (s+nz/2+zodd)*(size_t)nxy+(r+ny/2+yodd)*(size_t)nx+c+nx/2+xodd;
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
					idx1 = (size_t)s*nxy+(size_t)r*nx+c;
					idx2 = (s+nz/2+zodd)*(size_t)nxy+(r+ny/2+yodd)*(size_t)nx+c-nx/2-xodd;
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
					idx1 = (size_t)s*nxy+(size_t)r*nx+c;
					idx2 = (s+nz/2+zodd)*(size_t)nxy+(r-ny/2-yodd)*(size_t)nx+c-nx/2-xodd;
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
					idx1 = (size_t)s*nxy+(size_t)r*nx+c;
					idx2 = (s+nz/2+zodd)*(size_t)nxy+(r-ny/2-yodd)*(size_t)nx+c+nx/2+xodd;
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
					idx1 = (size_t)s*nxy+(size_t)r*nx+c;
					idx2 = (s+nz/2+zodd)*(size_t)nxy+(r+ny/2+yodd)*(size_t)nx+c;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}

			for( int s = nz-1; s >= (nz/2+zodd); --s ) {
				for ( int r = 0; r < ny/2; ++r ) {
					idx1 = (size_t)s*nxy+(size_t)r*nx+c;
					idx2 = (s-nz/2-zodd)*(size_t)nxy+(r+ny/2+yodd)*(size_t)nx+c;
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
					idx1 = (size_t)s*nxy+(size_t)r*nx+c;
					idx2 =(s+nz/2+zodd)*(size_t)nxy+(size_t)r*nx+c+nx/2+xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}

			for( int s = nz-1; s >= (nz/2+zodd); --s ) {
				for ( int c = 0; c < nx/2; ++c ) {
					idx1 = (size_t)s*nxy+(size_t)r*nx+c;
					idx2 = (s-nz/2-zodd)*(size_t)nxy+(size_t)r*nx+c+nx/2+xodd;
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
					idx1 = (size_t)s*nxy+(size_t)r*nx+c;
					idx2 = (size_t)s*nxy+(r+ny/2+yodd)*(size_t)nx+c+nx/2+xodd;
					tmp = rdata[idx1];
					rdata[idx1] = rdata[idx2];
					rdata[idx2] = tmp;
				}
			}

			for( int r = ny-1; r >= (ny/2+yodd); --r ) {
				for ( int c = 0; c < nx/2; ++c ) {
					idx1 = (size_t)s*nxy+(size_t)r*nx+c;
					idx2 = (size_t)s*nxy+(r-ny/2-yodd)*(size_t)nx+c+nx/2+xodd;
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

#ifdef EMAN2_USING_CUDA
	if (EMData::usecuda == 1 && image->getcudarwdata() && image->get_ndim() == 2) { // Because CUDA phase origin to center only works for 2D atm
		//cout << "CUDA tocorner " << image->getcudarwdata() << endl;
		emdata_phaseorigin_to_corner(image->getcudarwdata(), image->get_xsize(), image->get_ysize(), image->get_zsize());
		return;
	}
#endif // EMAN2_USING_CUDA

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
						idx = (size_t)s*nxy+r*nx+c;
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
						idx = (size_t)s*nxy+r*nx+c;
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
						idx = (size_t)s*nxy+r*nx+c;
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

#ifdef EMAN2_USING_CUDA
	if (EMData::usecuda == 1 && image->getcudarwdata() && image->get_ndim() == 2) { // Because CUDA phase origin to center only works for 2D atm
		//cout << "CUDA tocenter" << endl;
		emdata_phaseorigin_to_center(image->getcudarwdata(), image->get_xsize(), image->get_ysize(), image->get_zsize());
		return;
	}
#endif // EMAN2_USING_CUDA

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
						idx = (size_t)s*nxy+r*nx+c;
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
						idx = (size_t)s*nxy+r*nx+c;
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
						idx = (size_t)s*nxy+r*nx+c;
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

	if (image->get_ndim() != 2) {
		throw ImageDimensionException("This processor only supports 2D images.");
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

	int radius=params.set_default("radius",0);
	int nmaxseed=params.set_default("nmaxseed",0);
	
	float threshold=0.0;
	if (params.has_key("sigma") || !params.has_key("threshold")) threshold=(float)(image->get_attr("mean"))+(float)(image->get_attr("sigma"))*(float)params["sigma"];
	else threshold=params["threshold"];


	int nshells = params.set_default("nshells",0);
	int nshellsgauss = params.set_default("nshellsgauss",0);
	int verbose=params.set_default("verbose",0);

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	EMData *amask = new EMData();
	amask->set_size(nx, ny);

	float *dat = image->get_data();
	float *dat2 = amask->get_data();
	int i,j;
	size_t l = 0;

	if (verbose) printf("%f\t%f\t%f\t%d\n",(float)image->get_attr("mean"),(float)image->get_attr("sigma"),threshold,nmaxseed);

	// Seeds with the highest valued pixels
	if (nmaxseed>0) {
		EMData *peaks=image->process("mask.onlypeaks",Dict("npeaks",0));		// only find true peak values (pixels surrounded by lower values)
		vector<Pixel> maxs=peaks->calc_n_highest_locations(nmaxseed);
		delete peaks;
		
		for (vector<Pixel>::iterator i=maxs.begin(); i<maxs.end(); i++) {
			if ((*i).x==0 || (*i).y==0 ) continue;		// generally indicates a failed peak search, and regardless we don't really want edges
			amask->set_value_at((*i).x,(*i).y,0,1.0);
			if (verbose) printf("Seed at %d,%d,%d (%1.3f)\n",(*i).x,(*i).y,(*i).z,(*i).value);
		}
	}

	// Seeds with a circle
	if (radius>0) {
		// start with an initial circle
		l=0;
		for (j = -ny / 2; j < ny / 2; ++j) {
			for (i = -nx / 2; i < nx / 2; ++i,++l) {
				if ( abs(j) > radius || abs(i) > radius) continue;
				if ( (j * j + i * i) > (radius*radius) || dat[l] < threshold) continue;		// torn on the whole threshold issue here. Removing it prevents images from being totally masked out
				if ( (j * j + i * i) > (radius*radius) ) continue;
				dat2[l] = 1.0f;
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
		for (j=1; j<ny-1; ++j) {
			for (i=1; i<nx-1; ++i) {
				l=i+j*nx;
				if (dat2[l]) continue;
				if (dat[l]>threshold && (dat2[l-1]||dat2[l+1]||dat2[l+nx]||dat2[l-nx])) {
					dat2[l]=1.0;
					done=0;
				}
			}
		}
	}

	amask->update();

	if (verbose) printf("extending mask\n");
	if (nshells>0 || nshellsgauss>0) amask->process_inplace("mask.addshells.gauss", Dict("val1", nshells, "val2", nshellsgauss));

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

void CtfSimProcessor::process_inplace(EMData *image) {
	EMData *tmp=process(image);
	memcpy(image->get_data(),tmp->get_data(),(size_t)image->get_xsize()*image->get_ysize()*image->get_zsize()*sizeof(float));
	delete tmp;
	image->update();
	return;
}

EMData* CtfSimProcessor::process(const EMData * const image) {
	if (!image) {
		LOGWARN("NULL Image");
		return NULL;
	}

	EMData *fft;
	if (!image->is_complex()) fft=image->do_fft();
	else fft=image->copy();

	EMAN2Ctf ctf;
	ctf.defocus=params["defocus"];
	ctf.bfactor=params["bfactor"];
	ctf.ampcont=params.set_default("ampcont",10.0f);
	ctf.voltage=params.set_default("voltage",200.0f);
	ctf.cs=params.set_default("cs",2.0);
	ctf.apix=params.set_default("apix",image->get_attr_default("apix_x",1.0));
	ctf.dsbg=1.0/(ctf.apix*fft->get_ysize()*4.0);		//4x oversampling

	float noiseamp=params.set_default("noiseamp",0.0f);
	float noiseampwhite=params.set_default("noiseampwhite",0.0f);

	// compute and apply the CTF
	vector <float> ctfc = ctf.compute_1d(fft->get_ysize()*6,ctf.dsbg,ctf.CTF_AMP,NULL); // *6 goes to corner, remember you provide 2x the number of points you need

// 	printf("%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%d\n",ctf.defocus,ctf.bfactor,ctf.ampcont,ctf.dsbg,ctf.apix,fft->get_ysize());
// 	FILE *out=fopen("x.txt","w");
// 	for (int i=0; i<ctfc.size(); i++) fprintf(out,"%f\t%1.3g\n",0.25*i/(float)fft->get_ysize(),ctfc[i]);
// 	fclose(out);

	fft->apply_radial_func(0,0.25f/fft->get_ysize(),ctfc,1);

	// Add noise
	if (noiseamp!=0 || noiseampwhite!=0) {
		EMData *noise = new EMData(image->get_ysize(),image->get_ysize(),1);
		noise->process_inplace("testimage.noise.gauss");
		noise->do_fft_inplace();

		// White noise
		if (noiseampwhite!=0) {
			noise->mult((float)noiseampwhite*15.0f);		// The 15.0 is to roughly compensate for the stronger pink noise curve
			fft->add(*noise);
			noise->mult((float)1.0/(noiseampwhite*15.0f));
		}

		// Pink noise
		if (noiseamp!=0) {
			vector <float> pinkbg;
			pinkbg.resize(500);
			float nyimg=0.5f/ctf.apix;	// image nyquist
			// This pink curve came from a typical image in the GroEL 4A data set
			for (int i=0; i<500; i++) pinkbg[i]=noiseamp*(44.0f*exp(-5.0f*nyimg*i/250.0f)+10.0f*exp(-90.0f*nyimg*i/250.0f));		// compute curve to image Nyquist*2
			noise->apply_radial_func(0,.002f,pinkbg,1);		// Image nyquist is at 250 -> 0.5
			fft->add(*noise);
		}

	}

	EMData *ret=fft->do_ift();
	delete fft;

	return ret;
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

void IterMultiMaskProcessor::process_inplace(EMData * image)
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

	int num_shells = params.set_default("nshells",1);

	// there are other strategies which might allow us to avoid the extra copy, but this will have to do for now
	EMData *image1=image;
	EMData *image2=image->copy();
	if (nz == 1) {
		for (int i = 0; i < num_shells; i++) {
			for (int y = 1; y < ny - 1; y++) {
				for (int x = 1; x < nx - 1; x++) {
					if (image1->get_value_at(x,y)>0) continue;		// already part of a masked region

					// Note that this produces a directional bias in the case of ambiguous pixels
					// While this could be improved upon slightly, there can be truly ambiguous cases
					// and at least this method is deterministic
					if      (image1->get_value_at(x-1,y)>0) image2->set_value_at_fast(x,y,image1->get_value_at(x-1,y));
					else if (image1->get_value_at(x+1,y)>0) image2->set_value_at_fast(x,y,image1->get_value_at(x+1,y));
					else if (image1->get_value_at(x,y-1)>0) image2->set_value_at_fast(x,y,image1->get_value_at(x,y-1));
					else if (image1->get_value_at(x,y+1)>0) image2->set_value_at_fast(x,y,image1->get_value_at(x,y+1));

				}
			}
			memcpy(image1->get_data(),image2->get_data(),image1->get_size()*sizeof(float));
		}
	}
	else {
		for (int i = 0; i < num_shells; i++) {
			for (int z = 1; z < nz - 1; z++) {
				for (int y = 1; y < ny - 1; y++) {
					for (int x = 1; x < nx - 1; x++) {
						if (image1->get_value_at(x,y,z)>0) continue;		// already part of a masked region

						// Note that this produces a directional bias in the case of ambiguous pixels
						// While this could be improved upon slightly, there can be truly ambiguous cases
						// and at least this method is deterministic
						if      (image1->get_value_at(x-1,y,z)>0) image2->set_value_at_fast(x,y,z,image1->get_value_at(x-1,y,z));
						else if (image1->get_value_at(x+1,y,z)>0) image2->set_value_at_fast(x,y,z,image1->get_value_at(x+1,y,z));
						else if (image1->get_value_at(x,y-1,z)>0) image2->set_value_at_fast(x,y,z,image1->get_value_at(x,y-1,z));
						else if (image1->get_value_at(x,y+1,z)>0) image2->set_value_at_fast(x,y,z,image1->get_value_at(x,y+1,z));
						else if (image1->get_value_at(x,y,z-1)>0) image2->set_value_at_fast(x,y,z,image1->get_value_at(x,y,z-1));
						else if (image1->get_value_at(x,y,z+1)>0) image2->set_value_at_fast(x,y,z,image1->get_value_at(x,y,z+1));

					}
				}
			}
			memcpy(image1->get_data(),image2->get_data(),image1->get_size()*sizeof(float));
		}
	}

	delete image2;
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

	int num_shells = params.set_default("nshells",1);

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
				size_t cur_z = (size_t)z * nx * ny;

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

	size_t size = (size_t)nx * ny * nz;
	for (size_t i = 0; i < size; ++i) {
		if (d[i]) {
			d[i] = 1;
		}
		else {
			d[i] = 0;
		}
	}

	image->update();
}

void ToCenterProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	if ((float)image->get_attr("sigma")==0.0f) return;		// Can't center a constant valued image
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	// Preprocess a copy of the image to better isolate the "bulk" of the object
	float gmw=(nx/16)>5?nx/16:5;		// gaussian mask width
	EMData *image2=image->process("filter.highpass.gauss",Dict("cutoff_pixels",nx<50?nx/10:5));		// clear out large scale gradients
	image2->process_inplace("normalize.circlemean",Dict("radius",ny/2-4));
	image2->process_inplace("mask.gaussian",Dict("inner_radius",nx/2-gmw,"outer_radius",gmw/1.3));	// get rid of peripheral garbage
	image2->process_inplace("math.squared");		// exaggerate stronger density and includes strong negative density
	image2->process_inplace("filter.lowpass.gauss",Dict("cutoff_abs",0.05));						// get rid of peripheral garbage
	image2->process_inplace("normalize.circlemean",Dict("radius",ny/2-6));
	
	// We compute a histogram so we can decide on a good threshold value
	float hmin=(float)image2->get_attr("mean");
	float hmax=(float)image2->get_attr("mean")+(float)image2->get_attr("sigma")*4.0;
	vector <float> hist = image2->calc_hist( 100,hmin,hmax);
	double tot=0;
	int i;
	for (i=99; i>=0; i--) {
		tot+=hist[i];
		if (tot>nx*ny*nz/10) break;		// find a threshold encompassing a specific fraction of the total area/volume
	}
	float thr=(i*hmax+(99-i)*hmin)/99.0;		// this should now be a threshold encompasing ~1/10 of the area in the image
//	printf("mean %f   sigma %f   thr %f\n",(float)image2->get_attr("mean"),(float)image2->get_attr("sigma"),thr);
	
	// threshold so we are essentially centering the object silhouette
	image2->process_inplace("threshold.belowtozero",Dict("minval",thr));
//	image2->process_inplace("threshold.binary",Dict("value",thr));
//	image2->write_image("dbg1.hdf",-1);

	EMData *image3;
	if (nz==1) image3=image2->process("mask.auto2d",Dict("radius",nx/10,"threshold",thr*0.9,"nmaxseed",5));
	else image3=image2->process("mask.auto3d",Dict("radius",nx/10,"threshold",thr*0.9,"nmaxseed",5));

	// if the mask failed, we revert to the thresholded, but unmasked image
	if (nz==1 && (float)image3->get_attr("sigma")==0) {
		delete image3;
		image3=image2->process("math.linearpyramid");
		image3->add(9.0f);		// we comress the pyramid from .9-1
		image3->mult(0.9f);
		image3->process_inplace("mask.auto2d",Dict("threshold",0.5,"nmaxseed",5));	// should find seed points with a central bias
//		image3->write_image("dbg3.hdf",-1);
	}
		
	image3->process_inplace("threshold.binary",Dict("value",thr));

	if ((float)image3->get_attr("sigma")==0) delete image3;
	else {
		delete image2;
		image2=image3;
	}
	
//	image2->write_image("dbg2.hdf",-1);
	FloatPoint com = image2->calc_center_of_mass(0.5);
	delete image2;
	
	// actual centering is int translate only
	int dx = -(floor(com[0] + 0.5f) - nx / 2);
	int dy = -(floor(com[1] + 0.5f) - ny / 2);
	int dz = 0;
	if (nz > 1) {
		dz = -(floor(com[2] + 0.5f) - nz / 2);
	}
	if (abs(dx)>=nx-1 || abs(dy)>=ny-1 || abs(dz)>=nz) {
		printf("ERROR, center of mass outside image\n");
	}
	else {
		image->translate(dx, dy, dz);

		Transform t;
		t.set_trans((float)dx,(float)dy,(float)dz);

		if (nz > 1) {
			image->set_attr("xform.align3d",&t);
		} else {
			image->set_attr("xform.align2d",&t);
		}
	}

	
}


void ToMassCenterProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}

	int int_shift_only = params.set_default("int_shift_only",1);
	float threshold = params.set_default("threshold",0.0f);
//	int positive = params.set_default("positive",0);

	if ((float)image->get_attr("sigma")==0.0f) return;		// Can't center a constant valued image
	if (threshold>(float)image->get_attr("maximum")) {
		printf("Warning, centering threshold %1.2f, but image max %1.2f. Adjusting.",threshold,(float)image->get_attr("maximum"));
		threshold=(float)image->get_attr("mean")+(float)image->get_attr("sigma");
	}

	FloatPoint com = image->calc_center_of_mass(threshold);

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	if (int_shift_only) {
		int dx = -(floor(com[0] + 0.5f) - nx / 2);
		int dy = -(floor(com[1] + 0.5f) - ny / 2);
		int dz = 0;
		if (nz > 1) {
			dz = -(floor(com[2] + 0.5f) - nz / 2);
		}
		if (abs(dx)>=nx-1 || abs(dy)>=ny-1 || abs(dz)>=nz) {
			printf("ERROR, center of mass outside image\n");
		}
		else {
			image->translate(dx, dy, dz);

			Transform t;
			t.set_trans((float)dx,(float)dy,(float)dz);

			if (nz > 1) {
				image->set_attr("xform.align3d",&t);
			} else {
				image->set_attr("xform.align2d",&t);
			}
		}
	}
	else {
		float dx = -(com[0] - nx / 2);
		float dy = -(com[1] - ny / 2);
		float dz = 0;
		if (nz > 1) {
			dz = -(com[2] - nz / 2);
		}
		if (fabs(dx)>=nx-1 || fabs(dy)>=ny-2 || fabs(dz)>=nz) {
			printf("ERROR, center of mass outside image\n");
		}
		else {
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

void FSCFourierProcessor::process_inplace(EMData *image)
{
	EMData *tmp=process(image);
	size_t n = (size_t)image->get_xsize()*image->get_ysize()*image->get_zsize();
	memcpy(image->get_data(),tmp->get_data(),n*sizeof(float));
	image->update();
	delete tmp;
}

EMData *FSCFourierProcessor::process(EMData const *image)
{
	const char *fsp = params["fscfile"];
	float snrmult = params.set_default("snrmult",2.0f);
	float sscale = params.set_default("sscale",1.0f);
	float maxfreq = params.set_default("sscale",1.0f);

	XYData fsc;
	fsc.read_file(fsp);
	float nyquist=1.0/(2.0f*(float)image->get_attr("apix_x"));

	float lf=1.0f;
	int N=(int)fsc.get_size();
	int localav=0;				// once triggered, this uses a local average of 5 points instead of the point itself
	// While this could all be in one equation, the compiler will optimize it, and this is much more transparent
	for (int i=0; i<N; i++) {
		float s=i*nyquist/N;
		float f=fsc.get_y(i);
		float snr;
		if (s>=maxfreq && lf<f) f=lf;
		if (f<0 && i>2) localav=1;
		if (localav==1 && i>N-3) f=lf;
		else if (localav==1) f=(fsc.get_y(i-2)+fsc.get_y(i-1)+fsc.get_y(i)+fsc.get_y(i+1)+fsc.get_y(i+2))/5.0f;
		else if (localav==2) f=.00001;

		if (f>=1.0) f=.9999;
		if (f<0) { localav=2; f=.00001; }

		snr=snrmult*f/(1.0-f);		// if FSC==1, we just set it to 1000, which is large enough to make the Wiener filter effectively 1
		float wiener=snr*snr/(snr*snr+1);
		if (wiener<.001) wiener=.001;	// we don't want to go all the way to zero. We leave behind just just a touch to preserve potential phase info

		fsc.set_y(i,wiener);
		lf=f;
	}
	fsc.set_x(0,0);		// just to make sure we have values to the origin.
//	fsc.write_file("wiener.txt");
	FILE *out=fopen("wiener.txt","w");
	vector<float> wienerary(image->get_ysize());
	for (int i=0; i<image->get_ysize(); i++) {
		wienerary[i]=fsc.get_yatx(i*nyquist/image->get_ysize());
		fprintf(out,"%f\t%f\n",sscale*i*nyquist/image->get_ysize(),wienerary[i]);
	}
	fclose(out);

	EMData *fft=image->do_fft();
	fft->apply_radial_func(0,sscale*0.5/(float)image->get_ysize(),wienerary);

	EMData *ret=fft->do_ift();
	delete fft;

	return ret;
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

//	d2->apply_radial_func(0, 2.0f / Ctf::CTFOS, ctf, 0);	// this is troubling! only EMAN1 CTF used CTFOS, and then, not correctly...
	d2->apply_radial_func(0, 2.0f, ctf, 0);

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

void BadLineXYProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	if (image->get_zsize()>1) throw ImageDimensionException("Error: math.xybadline works only on 2D images");
	
	int xloc = params.set_default("xloc",-1);
	int yloc = params.set_default("yloc",-1);

	int nx=image->get_xsize();
	int ny=image->get_ysize();
	
	// each pixel is the average of its adjacent neighbors
	if (xloc>0 && xloc<nx-1) {
		for (int y=0; y<ny; y++) image->set_value_at(xloc,y,(image->get_value_at(xloc-1,y)+image->get_value_at(xloc+1,y))/2.0);
	}
	
	if (yloc>0 && yloc<ny-1) {
		for (int x=0; x<nx; x++) image->set_value_at(x,yloc,(image->get_value_at(x,yloc-1)+image->get_value_at(x,yloc+1))/2.0);
	}
}

void StripeXYProcessor::process_inplace(EMData * image)
{
	if (!image) {
		LOGWARN("NULL Image");
		return;
	}
	int xlen = params.set_default("xlen",10);
	int ylen = params.set_default("ylen",10);

	int nx=image->get_attr("nx");
	int ny=image->get_attr("ny");
	EMData *tmp=new EMData(nx,ny,1);

	// we do this in real-space, since the integration size is small, and we don't want Fourier edge effects
	// we use a 'moving window' to avoid summing the same points multiple times
	// start with y
	if (ylen>0) {
		for (int x=0; x<nx; x++) {
			float sum=0.0;
			float sumn=0.0;
			for (int y=0; y<(ylen<ny?ylen:ny); y++) { sum+=image->get_value_at(x,y); sumn+=1.0; } // we seed the process with a 1/2 stripe
			// now loop over each final pixel
			for (int y=0; y<ny; y++) {
				if (y+ylen<ny) { sum+=image->get_value_at(x,y+ylen); sumn+=1.0; }
				if (y-ylen-1>=0) { sum-=image->get_value_at(x,y-ylen-1); sumn-=1.0; }
				tmp->set_value_at_fast(x,y,sum/sumn);
			}
		}
		tmp->write_image("tmp.hdf",0);
		image->sub(*tmp);
	}

	// now x
	if (xlen>0) {
		for (int y=0; y<ny; y++) {
			float sum=0.0;
			float sumn=0.0;
			for (int x=0; x<(xlen<nx?xlen:nx); x++) { sum+=image->get_value_at(x,y); sumn+=1.0; } // we seed the process with a 1/2 stripe
			// now loop over each final pixel
			for (int x=0; x<nx; x++) {
				if (x+xlen<nx) { sum+=image->get_value_at(x+xlen,y); sumn+=1.0; }
				if (x-xlen-1>=0) { sum-=image->get_value_at(x-xlen-1,y); sumn-=1.0; }
				tmp->set_value_at_fast(x,y,sum/sumn);
			}
		}
		tmp->write_image("tmp.hdf",1);
		image->sub(*tmp);
	}

	delete tmp;
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
	maskblur->process_inplace("filter.lowpass.gauss", Dict("cutoff_pixels", radius));
//	maskblur->process_inplace("filter.highpass.tanh", Dict("highpass", -10.0f));
	maskblur->process_inplace("threshold.belowtozero", Dict("minval", 0.001f));
//	maskblur->process_inplace("threshold.belowtozero", Dict("minval", 0.001f));


	blur->process_inplace("threshold.belowtozero", Dict("minval", threshold));
	blur->process_inplace("filter.lowpass.gauss", Dict("cutoff_pixels", radius));
//	blur->process_inplace("filter.highpass.tanh", Dict("cutoff_abs", -10.0f));

	maskblur->div(*blur);
	image->mult(*maskblur);
//	maskblur->write_image("norm.mrc", 0, EMUtil::IMAGE_MRC);

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
				size_t index = (size_t)k * nx * ny + j * nx + i;
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

							size_t ii = x + y * nx + z * (size_t)xy;

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

	float xo = image->get_attr_default("origin_x",0.0);
	float yo = image->get_attr_default("origin_y",0.0);
	float zo = image->get_attr_default("origin_z",0.0);

	float xom = msk->get_attr_default("origin_x",0.0);
	float yom = msk->get_attr_default("origin_y",0.0);
	float zom = msk->get_attr_default("origin_z",0.0);

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
						idx = (size_t)k * nxy + j * nx + i;
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

void MatchSFProcessor::create_radial_func(vector < float >&rad,EMData *image) const {
	// The radial mask comes in with the existing radial image profile
	// The radial mask runs from 0 to the 1-D Nyquist (it leaves out the corners in Fourier space)

	EMData *to = params["to"];
	float apixto = to->get_attr("apix_x");
	bool bydot = params.set_default("bydot",false);
	float keephires = params.set_default("keephires",false);
	
	if (to->get_ysize() != image->get_ysize()) throw ImageDimensionException("Fatal Error: filter.matchto - image sizes must match");
	
	// This method is similar to that used by math.sub.optimal
	// It scales such that the average dot product at each resolution is matched,
	// but negative values (phase flips) are set to zero.
	if (bydot) {
		// Note that 'image' is guaranteed to be the FFT already
		EMData *tofft;
		if (to->is_complex()) tofft=to;
		else tofft=to->do_fft();
		
		float cornerscale;
		if (image->get_zsize()>1) cornerscale=sqrt(3.0);
		else cornerscale=sqrt(2.0);
		
		int ny=image->get_ysize();
		int ny2=(int)floor(image->get_ysize()*cornerscale/2);
		vector <float>norm(ny2+1);
		for (int i=0; i<ny2; i++) norm[i]=rad[i]=0;
		
		for (int y=-ny; y<ny; y++) {
			for (int x=0; x<ny; x++) {					// intentionally skipping the final pixel in x
				int r=int(Util::hypot_fast(x,y));
				if (r>ny2) continue;
				std::complex<float> v1=image->get_complex_at(x,y);
				std::complex<float> v2=tofft->get_complex_at(x,y);
				rad[r]+=(double)(v1.real()*v2.real()+v1.imag()*v2.imag());
				norm[r]+=(double)(v2.real()*v2.real()+v2.imag()*v2.imag());
			}
		}
		float radmean=0.0f;
		for (int i=0; i<ny2; i++) {
			rad[i]=(rad[i]/norm[i]<0)?0.0:rad[i]/norm[i];
			radmean+=rad[i];
		}
		radmean/=ny2;
		
		if (keephires) {
			for (int i=0; i<ny2; i++) {
				if (rad[i]<radmean/100.0) rad[i]=radmean/100.0;
			}
		}
			
		if (!to->is_complex()) delete tofft;
		return;
	}

	// This simply divides one radial intensity profile by the other
	if (to->is_complex()) {
		vector<float> rd=to->calc_radial_dist(rad.size(),0,1.0f,1);
		for (size_t i=0; i<rd.size(); ++i) {
			if (rad[i]>0) rad[i]=sqrt(rd[i]/rad[i]);
		}
	}
	else {
		EMData *tmp=to->do_fft();
		vector<float> rd=tmp->calc_radial_dist(rad.size(),0,1.0f,1);
		for (size_t i=0; i<rd.size(); ++i) {
			if (rad[i]>0) rad[i]=sqrt(rd[i]/rad[i]);
		}
		delete tmp;
	}
	
	// This makes sure that the filter never falls below 0.01 of the mean filter value so we keep something at all resolutions
	if (keephires) {
		float radmean=0.0f;
		for (int i=0; i<rad.size(); i++) radmean+=rad[i];
		radmean/=rad.size();
		for (int i=0; i<rad.size(); i++) {
			if (rad[i]<radmean/100.0) rad[i]=radmean/100.0;
		}
	}
		

}

void SetSFProcessor::create_radial_func(vector < float >&radial_mask,EMData *image) const {
	// The radial mask comes in with the existing radial image profile
	// The radial mask runs from 0 to the corner in Fourier space

	XYData *sf = params["strucfac"];
	if(params.has_key("apix")) {
		image->set_attr("apix_x", (float)params["apix"]);
		image->set_attr("apix_y", (float)params["apix"]);
		image->set_attr("apix_z", (float)params["apix"]);
	}

	float apix=image->get_attr("apix_x");
	int ny=image->get_ysize();
	int n = radial_mask.size();
	int nmax=(int)floor(sf->get_x(sf->get_size()-1)*apix*ny);		// This is the radius at which we have our last valid value from the curve
	if (nmax>n) nmax=n;

	if ((nmax)<3) throw InvalidParameterException("Insufficient structure factor data for SetSFProcessor to be meaningful");

	int i;
	for (i=0; i<nmax; i++) {
//		if (radial_mask[i]>0)
//		{
//			radial_mask[i]= sqrt(n*sf->get_yatx(i/(apix*2.0f*n),false)/radial_mask[i]);
//		}
		if (radial_mask[i]>0) {
			radial_mask[i]= sqrt((ny*ny*ny)*sf->get_yatx(i/(apix*ny))/radial_mask[i]);
		}
		else if (i>0) radial_mask[i]=radial_mask[i-1];	// For points where the existing power spectrum was 0
	}

	// Continue to use a fixed factor after we run out of 'sf' values

	while (i<n) {
		radial_mask[i]=radial_mask[nmax-1];
		i++;
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
	size_t smn = 0;
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

	float smask = (float) (sma / smn);
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
			size_t k2 = (size_t)k * nxy;
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
	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			int j2 = j * nx + i;
			int k0 = 0;
			for (int k = 0; k < nz; ++k) {
				idx = j2 + (size_t)k * nxy;
				if (dat2[idx]) {
					k0 = k;
					break;
				}
			}

			if (k0 != nz) {
				int k1 = nz - 1;
				for (int k = nz - 1; k >= 0; --k) {
					idx = j2 + (size_t)k * nxy;
					if (dat2[idx]) {
						k1 = k;
						break;
					}
				}

				for (int k = k0 + 1; k < k1; ++k) {
					idx = j2 + (size_t)k * nxy;
					dat2[idx] = 1.0f;
				}
			}
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < nz; ++j) {
			size_t j2 = (size_t)j * nxy + i;
			int k0 = 0;
			for (int k = 0; k < ny; ++k) {
				idx = (size_t)k * nx + j2;
				if (dat2[idx]) {
					k0 = k;
					break;
				}
			}

			if (k0 != ny) {
				int k1 = ny - 1;
				for (int k = ny - 1; k >= 0; --k) {
					idx = (size_t)k * nx + j2;
					if (dat2[idx]) {
						k1 = k;
						break;
					}
				}

				for (int k = k0 + 1; k < k1; ++k) {
					idx = (size_t)k * nx + j2;
					dat2[idx] = 1.0f;
				}
			}
		}
	}

	for (int i = 0; i < ny; ++i) {
		for (int j = 0; j < nz; ++j) {
			size_t j2 = i * nx + (size_t)j * nxy;
			int k0 = 0;
			for (int k = 0; k < nx; ++k) {
				if (dat2[k + j2]) {
					k0 = k;
					break;
				}
			}
			if (k0 != nx) {
				int k1 = nx - 1;
				for (int k = nx - 1; k >= 0; --k) {
					if (dat2[k + j2]) {
						k1 = k;
						break;
					}
				}

				for (int k = k0 + 1; k < k1; ++k) {
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
			size_t l = j * nx + (size_t)k * nxy + 1;
			for (int i = 1; i < nx - 1; ++i, ++l) {
				if (dat2[l - 1] == 1.0f || dat2[l + 1] == 1.0f ||
					dat2[l - nx] == 1.0f || dat2[l + nx] == 1.0f ||
					dat2[l - nxy] == 1.0f || dat2[l + nxy] == 1.0f) {
					dat2[l] = 2.0f;
				}
			}
		}
	}

	size_t size = (size_t)nx * ny * nz;
	for (size_t i = 0; i < size; ++i) {
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

// Originally this was done recursively, but there were issues with exceeding the maximum recursion depth,
// so switched to a vector-style algorithm instead
void AutoMaskDustProcessor::process_inplace(EMData * imagein)
{
	if (!imagein) {
		LOGWARN("NULL Image");
		return;
	}
	image=imagein;

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	int verbose=params.set_default("verbose",0);
	unsigned int voxels=params.set_default("voxels",27);
	float threshold=params.set_default("threshold",1.5);

	mask = new EMData();
	mask->set_size(nx, ny, nz);
	mask->to_one();

	for (int zz = 0; zz < nz; zz++) {
		for (int yy = 0; yy < ny; yy++) {
			for (int xx = 0; xx < nx; xx++) {
				if (image->get_value_at(xx,yy,zz)>threshold && mask->get_value_at(xx,yy,zz)==1.0) {
					vector<Vec3i> pvec;
					pvec.push_back(Vec3i(xx,yy,zz));
					for (uint i=0; i<pvec.size(); i++) {
						// Duplicates will occur the way the algorithm is constructed, so we eliminate them as we encounter them
						if (mask->sget_value_at(pvec[i])==0.0f) {
							pvec.erase(pvec.begin()+i);
							i--;
							continue;
						}

						// mask out the points in the volume
						mask->set_value_at(pvec[i],0.0f);

						int x=pvec[i][0];
						int y=pvec[i][1];
						int z=pvec[i][2];
						// Any neighboring values above threshold we haven't already set get added to the list
						if (image->sget_value_at(x-1,y,z)>threshold && mask->get_value_at(x-1,y,z)==1.0) pvec.push_back(Vec3i(x-1,y,z));
						if (image->sget_value_at(x+1,y,z)>threshold && mask->get_value_at(x+1,y,z)==1.0) pvec.push_back(Vec3i(x+1,y,z));
						if (image->sget_value_at(x,y-1,z)>threshold && mask->get_value_at(x,y-1,z)==1.0) pvec.push_back(Vec3i(x,y-1,z));
						if (image->sget_value_at(x,y+1,z)>threshold && mask->get_value_at(x,y+1,z)==1.0) pvec.push_back(Vec3i(x,y+1,z));
						if (image->sget_value_at(x,y,z-1)>threshold && mask->get_value_at(x,y,z-1)==1.0) pvec.push_back(Vec3i(x,y,z-1));
						if (image->sget_value_at(x,y,z+1)>threshold && mask->get_value_at(x,y,z+1)==1.0) pvec.push_back(Vec3i(x,y,z+1));
					}

					// If the blob is too big, then we don't mask it out after all, but we set the value
					// to 2.0 so we know the voxels have already been examined, and don't check them again
					if (pvec.size()>voxels) {
						if (verbose) printf("%d\t%d\t%d\tvoxels: %d\n",xx,yy,zz,(int)pvec.size());
						for (uint i=0; i<pvec.size(); i++) mask->set_value_at(pvec[i],2.0);
					}
				}
			}
		}
	}

	mask->process_inplace("threshold.binary",Dict("value",0.5));

	// Now we expand the mask by 1 pixel and blur the edge
	mask->mult(-1.0f);
	mask->add(1.0f);
	mask->process_inplace("mask.addshells",Dict("nshells",2));		// expand by 1 shell
	mask->process_inplace("filter.lowpass.gauss",Dict("cutoff_abs",0.25f));
	mask->mult(-1.0f);
	mask->add(1.0f);
	mask->update();

	// apply the mask
	image->mult(*mask);
	if (verbose>1) mask->write_image("mask.hdf", 0, EMUtil::IMAGE_HDF);

	delete mask;
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

	float threshold=0.0;
	if (params.has_key("sigma")) threshold=(float)(image->get_attr("mean"))+(float)(image->get_attr("sigma"))*(float)params["sigma"];
	else threshold=params["threshold"];

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

	if (verbose) printf("expanding mask\n");
	amask->process_inplace("mask.addshells.gauss", Dict("val1", (int)(nshells+nshellsgauss/2),"val2",0));
	if (verbose) printf("filtering mask\n");
	amask->process_inplace("filter.lowpass.gauss", Dict("cutoff_abs", (float)(1.0f/(float)nshellsgauss)));
	amask->process_inplace("threshold.belowtozero", Dict("minval",(float)0.002));	// this makes the value exactly 0 beyond ~2.5 sigma

	bool return_mask = params.set_default("return_mask",false);
	if (return_mask) {
		// Yes there is probably a much more efficient way of getting the mask itself, but I am only providing a stop gap at the moment.
		dat = image->get_data();
		dat2 = amask->get_data();
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
	size_t size = (size_t)nx * ny * nz;

	// TODO: THIS IS EXTREMELY INEFFICIENT
	if (nz != 1) {
		for (int l = 1; l <= (int) val1+val2; ++l) {
			for (size_t i=0; i<size; i++) d2[i]=d[i];
			for (int k = 1; k < nz - 1; ++k) {
				for (int j = 1; j < ny - 1; ++j) {
					for (int i = 1; i < nx - 1; ++i) {
						size_t t = i + j*nx+(size_t)k*nx*ny;
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
	for (size_t i = 0; i < size; ++i) if (d[i]) d[i]=vec[(int)d[i]];

	image->update();
	delete image2;
}

EMData* DirectionalSumProcessor::process(const EMData* const image ) {
	string dir = params.set_default("axis", "");
	if ( dir == "" || ( dir != "x" && dir != "y" && dir != "z" ) )
		throw InvalidParameterException("The direction parameter must be either x, y, or z");

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	EMData* ret = new EMData;
	
	if (nz==1) {
		if (dir=="x") {
			ret->set_size(ny,1,1);
			for (int y=0; y<ny; y++) {
				double sm=0;
				for (int x=0; x<nx; x++) sm+=image->get_value_at(x,y);
				ret->set_value_at(y,0,0,sm/nx);
			}
			return ret;
		}
		else if (dir=="y") {
			ret->set_size(nx,1,1);
			for (int x=0; x<nx; x++) {
				double sm=0;
				for (int y=0; y<ny; y++) sm+=image->get_value_at(x,y);
				ret->set_value_at(x,0,0,sm/ny);
			}
			return ret;
		}
		else throw InvalidParameterException("The direction parameter must be either x, y for 2-D images");
	}
	
	int a0 = params.set_default("first", 0);
	int a1 = params.set_default("last", -1);

	// compress one of the dimensions
	if ( dir == "x" ) {
		ret->set_size(nz,ny);

		// bounds checks
		if (a0<0) a0+=nx;
		if (a1<0) a1+=nx;
		if (a0<0) a0=0;
		if (a1<0) a1=0;
		if (a0>=nx) a0=nx-1;
		if (a1>=nx) a1=nx-1;

		for (int y=0; y<ny; y++) {
			for (int z=0; z<nz; z++) {
				double sum=0.0;
				for (int x=a0; x<=a1; x++) sum+=image->get_value_at(x,y,z);
				ret->set_value_at(z,y,(float)sum);
			}
		}
	}
	else if ( dir == "y" ) {
		ret->set_size(nx,nz);

		// bounds checks
		if (a0<0) a0+=ny;
		if (a1<0) a1+=ny;
		if (a0<0) a0=0;
		if (a1<0) a1=0;
		if (a0>=ny) a0=ny-1;
		if (a1>=ny) a1=ny-1;

		for (int x=0; x<nx; x++) {
			for (int z=0; z<nz; z++) {
				double sum=0.0;
				for (int y=a0; y<=a1; y++) sum+=image->get_value_at(x,y,z);
				ret->set_value_at(x,z,(float)sum);
			}
		}
	}
	else if ( dir == "z" ) {
		ret->set_size(nx,ny);

		// bounds checks
		if (a0<0) a0+=nz;
		if (a1<0) a1+=nz;
		if (a0<0) a0=0;
		if (a1<0) a1=0;
		if (a0>=nz) a0=nz-1;
		if (a1>=nz) a1=nz-1;

		for (int y=0; y<ny; y++) {
			for (int x=0; x<nx; x++) {
				double sum=0.0;
				for (int z=a0; z<=a1; z++) sum+=image->get_value_at(x,y,z);
				ret->set_value_at(x,y,(float)sum);
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
				size_t idx1 = twox + y*nx+(size_t)z*nxy;
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
				size_t idx1 = twox + y*nx+(size_t)z*nxy;
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
				size_t idx1 = twox + y*nx+(size_t)z*nxy;
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
	size_t n = (size_t)image->get_xsize()*image->get_ysize()*image->get_zsize();

	for(size_t i = 0; i < n; ++i) {
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
	int fill = 1;
	if(params.has_key("fill")) {
		fill = (int)params["fill"];
	}

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
	size_t size = (size_t)nx*ny*nz;
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
	size_t size = (size_t)nx*ny*nz;
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

	if ((int)params["dir"]==1) dir=(gsl_wavelet_direction)1;	// defined as 'forward', but name is broken on mac
	else dir=(gsl_wavelet_direction)-1;	// backward

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
	} else if (axis == "y" || axis == "Y") {
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
	} else if (axis == "z" || axis == "Z") {
		if(1-z_start) {
			int nhalf = nz/2;
			float *tmp = new float[nxy];
			for(int iz = 0;iz<nhalf;iz++){
				memcpy(tmp,&data[iz*nxy],nxy*sizeof(float));
				memcpy(&data[iz*nxy],&data[(nz-iz-1)*nxy],nxy*sizeof(float));
				memcpy(&data[(nz-iz-1)*nxy],tmp,nxy*sizeof(float));
			}
			delete[] tmp;
		} else {
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
	int N   = ny;

	int zerocorners = params.set_default("zerocorners",0);

	
	const float * const src_data = image->get_const_data();
	float *des_data = (float *) EMUtil::em_malloc(nx*ny*nz* sizeof(float));

	if ((nz == 1)&&(image -> is_real()))  {
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
	if ((nz == 1)&&(image -> is_complex())&&(nx%2==0)&&((2*(nx-ny)-3)*(2*(nx-ny)-3)==1)&&(zerocorners==0) )  { 
 	  //printf("Hello 2-d complex  TransformProcessor \n");
	  // make sure there was a realImage.process('xform.phaseorigin.tocorner')
	 //           before transformation to Fourier Space
// This rotates a complex image that is a FT of a real space image: it has Friedel symmetries
// An arbitrary complex image, F, can be decomposed into two Friedel images 
//         G(k) = (F(k) + F*(-k))/2 ,  H(k) = (-i(F(k) - F*(k)))/2;
//          via   F(k) = G(k) + i H(k); notice G and H are Friedel symmetric (not necessarily real)
//	           But g,h are real, using f=g+ih.
// 		First make sure that image has proper size; 
//         if 2N is side of real image, then sizes of FFT are (2N+2,2N)
//         if 2N+1 is side of real image, then sizes of FFT are (2N+2,2N+1)
//         so we need nx =ny+2, and ny  even  
//	          or  nx = ny+1  and ny odd
//         so nx even, and 2 *(nx -ny) -3=  +- 1; So  abs(2*(nx-ny)-3) == 1
		float theta =  t.get_rotation("eman").get("phi"); theta=theta*pi/180;
		Vec3f transNow= t.get_trans();
		float xshift= transNow[0]; float yshift= transNow[1];
		float tempR; float tempI;float tempW;
//		int kNy= ny; //size of the real space image
//		int kNx= nx/2; //
		Vec2f offset(nx/2,ny/2); 
		float Mid =(N+1.0)/2.0;  // Check this
		for (int kyN = 0; kyN < ny; kyN++) {
			int kyNew = kyN;
			if (kyN>=nx/2) kyNew=kyN-ny;  // Step 0     Unalias
			float kxOldPre=  - sin(theta)* kyNew;
			float kyOldPre=    cos(theta)* kyNew;
    		        float phaseY = -2*pi*kyNew*yshift/ny;

			for (int kxN = 0; kxN < (nx/2); kxN++) {
				int kxNew=kxN; 
				if (kxN >= nx/2) kxNew=kxN-ny;//  Step 0   Unalias
				//  Step 1, Do rotation and find 4 nn
				float kxOld=  kxOldPre + cos(theta)* kxNew ;
				float kyOld=  kyOldPre + sin(theta)* kxNew ;
				//
				int kxLower= floor(kxOld); int kxUpper= kxLower+1;    // All these values
				int kyLower= floor(kyOld); int kyUpper= kyLower+1;    // may be neg
				float dkxLower= (kxUpper-kxOld);    float dkxUpper= (kxOld-kxLower);
				float dkyLower= (kyUpper-kyOld);    float dkyUpper= (kyOld-kyLower);
//
				int kxL= kxLower; int kyL=kyLower; 
				float dataLL_R= 0; float dataLL_I=0; int flag=1;
				if ((abs(kxL)<Mid) && (abs(kyL)<Mid)) { //   Step 2 Make sure to be in First BZ
				    kxL = (N+kxL)%N;  kyL = (N+kyL)%N;
				    if (kxL> N/2){ kxL=(N-kxL)%N; kyL=(N-kyL)%N ;flag=-1;} // Step 3: if nec, use Friedel paired
				    dataLL_R=      src_data[2*kxL+   kyL*nx]; //  image -> get_value_at(2*kxL,kyL);
				    dataLL_I= flag*src_data[2*kxL+1+ kyL*nx]; //  image -> get_value_at(2*kxL+1,kyL);
				} 

			        kxL=kxLower; int kyU=kyUpper; 
			        float dataLU_R= 0; float dataLU_I=0; flag=1;
				if ((abs(kxL)<Mid) && (abs(kyU)<Mid)){ //       Step 2 Make sure to be in First BZ
				    kxL = (N+kxL)%N;  kyU = (N+kyU)%N;
				    if (kxL> N/2){ kxL=(N-kxL)%N; kyU=(N-kyU)%N;flag=-1;} // Step 3
				    dataLU_R=	    src_data[2*kxL+   kyU*nx];//  image -> get_value_at(2*kxL,kyU);
				    dataLU_I=  flag*src_data[2*kxL+1+ kyU*nx];// flag*image -> get_value_at(2*kxL+1,kyU);
			        }

				int kxU= kxUpper; kyL=kyLower; 
				float dataUL_R= 0; float dataUL_I=0; flag=1;
				if ((abs(kxU)<Mid) && (abs(kyL)<Mid)) {   //       Step 2
				    kxU = (N+kxU)%N; kyL = (N+kyL)%N;
				    if (kxU> N/2) { kxU=(N-kxU)%N; kyL=(N-kyL)%N;flag=-1;} // Step 3
				    dataUL_R=	    src_data[2*kxU   +   kyL*nx]; //  image -> get_value_at(2*kxU,kyL);
				    dataUL_I=  flag*src_data[2*kxU+1 +   kyL*nx]; // flag*image -> get_value_at(2*kxU+1,kyL);
				}

			      kxU= kxUpper; kyU=kyUpper; 
			      float dataUU_R= 0; float dataUU_I=0; flag=1;
			      if ((abs(kxU)<Mid) & (abs(kyU)<Mid)){  //       Step 2
				    kxU = (N+kxU)%N; kyU = (N+kyU)%N;
				    if (kxU> N/2) { kxU=(N-kxU)%N; kyU=(N-kyU)%N;flag=-1;} // Step 3
				    dataUU_R=	    src_data[2*kxU   +   kyU*nx]; //   image -> get_value_at(2*kxU,kyU);
				    dataUU_I=  flag*src_data[2*kxU+1 +   kyU*nx]; //flag*image -> get_value_at(2*kxU+1,kyU);
			      }
			      //            Step 4    Assign Weights
			      float WLL = dkxLower*dkyLower  ;
			      float WLU = dkxLower*dkyUpper  ;// make more intricated weightings here 
			      float WUL = dkxUpper*dkyLower  ;// WLL(dkxLower,dkyLower)
			      float WUU = dkxUpper*dkyUpper  ;//  etc
			      tempW = WLL  +  WLU + WUL +  WUU ;

			      //            Step 5    Assign Real, then Imaginar Values
			      tempR = WLL*dataLL_R  +   WLU*dataLU_R + WUL* dataUL_R +   WUU * dataUU_R ;
			      tempI = WLL*dataLL_I  +   WLU*dataLU_I + WUL* dataUL_I +   WUU * dataUU_I ;
			      //
			      float phase = phaseY -2*pi*kxNew*xshift/ny;
			      float tempRb=tempR*cos(phase) - tempI*sin(phase);
			      float tempIb=tempR*sin(phase) + tempI*cos(phase);
			      //
			      des_data[2*kxN   + nx* kyN] = tempRb/tempW;
			      des_data[2*kxN+1 + nx* kyN] = tempIb/tempW;
			      //printf(" kxNew = %d, kyNew = %d,kxOld = %3.2f, kyOld = %3.2f,  xl = %d,xU = %d,yl = %d,yu = %d, tempR = %3.2f, tempI=%3.2f,  \n", 
				//      kxNew,kyNew, kxOld, kyOld, kxLower,kxUpper,kyLower,kyUpper, tempR, tempI);
			}
		}
	}

	if ((nz == 1)&&(image -> is_complex())&&(nx%2==0)&&((2*(nx-ny)-3)*(2*(nx-ny)-3)==1)&&(zerocorners==1) )  { 
 	  //printf("Hello 2-d complex  TransformProcessor \n");
	  // make sure there was a realImage.process('xform.phaseorigin.tocorner')
	 //           before transformation to Fourier Space
// This rotates a complex image that is a FT of a real space image: it has Friedel symmetries
// An arbitrary complex image, F, can be decomposed into two Friedel images 
//         G(k) = (F(k) + F*(-k))/2 ,  H(k) = (-i(F(k) - F*(k)))/2;
//          via   F(k) = G(k) + i H(k); notice G and H are Friedel symmetric (not necessarily real)
//	           But g,h are real, using f=g+ih.
// 		First make sure that image has proper size; 
//         if 2N is side of real image, then sizes of FFT are (2N+2,2N)
//         if 2N+1 is side of real image, then sizes of FFT are (2N+2,2N+1)
//         so we need nx =ny+2, and ny  even  
//	          or  nx = ny+1  and ny odd
//         so nx even, and 2 *(nx -ny) -3=  +- 1; So  abs(2*(nx-ny)-3) == 1
		float theta =  t.get_rotation("eman").get("phi"); theta=theta*pi/180;
		float Ctheta= cos(theta);
		float Stheta= sin(theta);
		Vec3f transNow= t.get_trans();
		float xshift= transNow[0]; float yshift= transNow[1];
		float tempR; float tempI; // float tempW;
		float Mid =(N+1.0)/2.0;  // Check this
//		int kNy= ny; //size of the real space image
//		int kNx= nx/2; //
		Vec2f offset(nx/2,ny/2); 
		float phaseConstx  = -2*pi*xshift/ny ;
		float k1= cos(phaseConstx); float k2= sin(phaseConstx);
		float k3= 1.0/k1; float k4= k2/k1; // that is 1/cos and tan()
		
		for (int kyN = 0; kyN < ny; kyN++) {
			int kyNew = kyN;
			if (kyN>=nx/2) kyNew=kyN-ny;  // Step 0     Unalias
			float kxOld=  - Stheta* kyNew - Ctheta;
			float kyOld=    Ctheta* kyNew - Stheta;
    		        float phase = -2*pi*kyNew*yshift/ny - phaseConstx ;
			float Cphase = cos(phase);
			float Sphase = sin(phase);
			int kx,ky;
			int IndexOut=  -2+ nx* kyN;
			for (int kxN = 0; kxN < (nx/2); kxN++) {
				IndexOut += 2;  
				// int kxNew=kxN; 
				// if (kxN >= nx/2) kxNew=kxN-ny;//  Step 0   Unalias, but this won't happen
				//  Step 1, Do rotation and find 4 nn
				kxOld +=  Ctheta ;// was using kxNew instead of kxN
				kyOld +=  Stheta ;
				phase += phaseConstx;
			        Cphase = Cphase*k1 -Sphase*k2; //update using trig addition; this is   cos = cos cos  -sin sin
			        Sphase = Sphase*k3+ Cphase*k4;  //   and   sin = sin  (1/ cos) + cos * tan;

				if ((abs(kxOld)>=Mid) || (abs(kyOld)>=Mid)) { // out of bounds
         			      des_data[IndexOut] = 0;
	        		      des_data[IndexOut+1] = 0;
				   continue;}
				   
				//
				int kxLower= floor(kxOld);     // All these values
				int kyLower= floor(kyOld);    // may be neg
				//float dkxLower= (kxLower+1-kxOld);    
				float dkxUpper= (kxOld-kxLower);
				//float dkyLower= (kyLower+1-kyOld);    
				float dkyUpper= (kyOld-kyLower);
//
				kx= kxLower; ky =kyLower; 
				int flagLL=1;
				    if (kx<0) kx += N;  
				    if (ky<0) ky += N;
				    if (kx> N/2){ kx=(N-kx) ; ky=(N-ky)  ;flagLL=-1;} // Step 3: if nec, use Friedel paired
				    int kLL =2*kx+ky*nx; 

			        kx = kxLower; ky = kyLower+1; 
				int flagLU =1;
				    if (kx<0) kx += N;  
				    if (ky<0) ky += N;
				    if (kx> N/2){ kx=(N-kx) ; ky=(N-ky)  ;flagLU=-1;} // Step 3: if nec, use Friedel paired
				    int kLU =2*kx+ky*nx; 

				kx = kxLower+1; ky=kyLower; 
				int flagUL=1;
				    if (kx<0) kx += N;  
				    if (ky<0) ky += N;
				    if (kx> N/2){ kx=(N-kx) ; ky=(N-ky)  ;flagUL=-1;} // Step 3: if nec, use Friedel paired
				    int kUL =2*kx+ky*nx; 

			        kx = kxLower+1; ky = kyLower+1; 
			        int flagUU =1;
				    if (kx<0) kx += N;  
				    if (ky<0) ky += N;
				    if (kx> N/2){ kx=(N-kx) ; ky=(N-ky)  ;flagUU=-1;} // Step 3: if nec, use Friedel paired
				    int kUU =2*kx+ky*nx; 
				    
			      //            Step 4    Assign Weights
// 			      float WLL = dkxLower*dkyLower  ;
// 			      float WLU = dkxLower*dkyUpper  ;// make more intricated weightings here 
// 			      float WUL = dkxUpper*dkyLower  ;// WLL(dkxLower,dkyLower)
// 			      float WUU = dkxUpper*dkyUpper  ;//  etc
// 			      tempW = WLL  +  WLU + WUL +  WUU ;

			      //            Step 5    Assign Real, then Imaginary Values
			      tempR = Util::bilinear_interpolate(src_data[kLL],src_data[kUL], src_data[kLU], src_data[kUU],dkxUpper,dkyUpper);//
			               //WLL*dataLL_R  +   WLU*dataLU_R + WUL* dataUL_R +   WUU * dataUU_R ;
			      tempI = Util::bilinear_interpolate(flagLL* src_data[kLL+1],flagUL*src_data[kUL+1], 
					                   flagLU*src_data[kLU+1], flagUU*src_data[kUU+1],dkxUpper,dkyUpper);
			               //WLL*dataLL_I  +   WLU*dataLU_I + WUL* dataUL_I +   WUU * dataUU_I ;
			      //
			      
			      float tempRb=tempR*Cphase - tempI*Sphase;
			      float tempIb=tempR*Sphase + tempI*Cphase;
			      //
			      des_data[IndexOut]   = tempRb;
			      des_data[IndexOut+1] = tempIb;
			      //printf(" kxNew = %d, kyNew = %d,kxOld = %3.2f, kyOld = %3.2f,  xl = %d,xU = %d,yl = %d,yu = %d, tempR = %3.2f, tempI=%3.2f,  \n", 
				//      kxNew,kyNew, kxOld, kyOld, kxLower,kxUpper,kyLower,kyUpper, tempR, tempI);
			}
		}
	}

	if ((nz > 1)&&(image -> is_complex())&&(zerocorners<2)) {
		//printf("Hello 3-d complex  TransformProcessor \n");
		float phi =  t.get_rotation("eman").get("phi"); phi=pi*phi/180;
		float alt =  t.get_rotation("eman").get("alt"); alt=pi*alt/180;
		float az  =  t.get_rotation("eman").get("az");   az=pi*az /180;
		Vec3f transNow= t.get_trans();
		float xshift= transNow[0]; float yshift= transNow[1];float zshift= transNow[2];

		float phaseConstx  = -2*pi*xshift/ny ;
		float k1= cos(phaseConstx); float k2= sin(phaseConstx);
		float k3= 1.0/k1; float k4= k2/k1; // that is 1/cos and tan()
		
		float MatXX = (cos(az)*cos(phi) - sin(az)*cos(alt)*sin(phi) );
		float MatXY = (- cos(az)*sin(phi) - sin(az)*cos(alt)*cos(phi) ) ;
		float MatXZ = sin(az)*sin(alt) ;
		float MatYX = (sin(az)*cos(phi) + cos(az)*cos(alt)*sin(phi) );
		float MatYY = (- sin(az)*sin(phi) + cos(az)*cos(alt)*cos(phi) )  ;
		float MatYZ = - cos(az)*sin(alt) ;
		float MatZX = sin(alt)*sin(phi);
		float MatZY = sin(alt)*cos(phi);
		float MatZZ = cos(alt)  ;
		float tempR; float tempI; float tempW;
		float Mid =(N+1.0)/2.0;  // Check this
		int nxny = nx*ny;
		
		for (int kzN = 0; kzN < ny; kzN++) { 
		    int kzNew=kzN; 
		    if (kzN >= nx/2) kzNew=kzN-N; //            Step 0 Unalias new coords; 
		    for (int kyN = 0; kyN < ny; kyN++) {//       moves them to lesser mag   
			int kyNew=kyN; 
			if (kyN>=nx/2) kyNew=kyN-ny;  //        Step 0   Unalias
		        float kxPre =  MatXY * kyNew +  MatXZ *kzNew;
	                float kyPre =  MatYY * kyNew +  MatYZ*kzNew;
			float kzPre =  MatZY * kyNew +  MatZZ*kzNew;
    		        float phase = -2*pi*kzNew*zshift/ny-2*pi*kyNew*yshift/ny - phaseConstx ;
			float Cphase = cos(phase);
			float Sphase = sin(phase);
			
			float OutBounds2= (2*Mid*Mid- (kyNew*kyNew+kzNew*kzNew)) ;
			int kxNewMax= nx/2;
			if (OutBounds2< 0) kxNewMax=0;
			else if (OutBounds2<(nx*nx/4)) kxNewMax=sqrt(OutBounds2);
			for (int kxN = kxNewMax; kxN < nx/2 ; kxN++ ) {
       			      des_data[2*kxN    + nx* kyN +nxny*kzN] = 0;
        		      des_data[2*kxN    + nx* kyN +nxny*kzN+1] = 0;
			}
			
			
			for (int kxN = 0; kxN < kxNewMax; kxN++ ) {
			     Cphase = Cphase*k1 -Sphase*k2; //update using trig addition; this is   cos = cos cos  -sin sin
			     Sphase = Sphase*k3+ Cphase*k4;  //   and   sin = sin  (1/ cos) + cos * tan;
			    //  Step 1: Do inverse Rotation to get former values, and alias   Step1
			    float kxOld=  MatXX * kxN + kxPre;
			    float kyOld=  MatYX * kxN + kyPre;
			    float kzOld=  MatZX * kxN + kzPre;
			    //
				if ((abs(kxOld)>=Mid) || (abs(kyOld)>=Mid) || (abs(kzOld)>=Mid) ) { // out of bounds
         			      des_data[2*kxN    + nx* kyN +nxny*kzN] = 0;
	        		      des_data[2*kxN    + nx* kyN +nxny*kzN+1] = 0;
				   continue;}
				//
			    int kxLower= floor(kxOld); int kxUpper= kxLower+1;
			    int kyLower= floor(kyOld); int kyUpper= kyLower+1;
			    int kzLower= floor(kzOld); int kzUpper= kzLower+1;
			    //
			    float dkxLower= (kxUpper-kxOld); float dkxUpper= (kxOld-kxLower);
			    float dkyLower= (kyUpper-kyOld); float dkyUpper= (kyOld-kyLower);
			    float dkzLower= (kzUpper-kzOld); float dkzUpper= (kzOld-kzLower);
			    //      LLL  1
			    int kxL= kxLower; int kyL=kyLower; int kzL=kzLower; 
			    float dataLLL_R= 0; float dataLLL_I=0; int flag=1;
			    if ( (abs(kxL)<Mid) && (abs(kyL)<Mid) && (abs(kzL)<Mid) ) { //   Step 2 Make sure to be in First BZ
				kxL = (N+kxL)%N;  kyL = (N+kyL)%N; kzL = (N+kzL)%N;
				if (kxL> N/2){kxL=(N-kxL)%N; kyL=(N-kyL)%N ; kzL=(N-kzL)%N ;flag=-1;} // Step 3: use Friedel paired
				dataLLL_R=     src_data[ 2*kxL   + nx*kyL+ nxny*kzL ]; //get_value_at(2*kxL,kyL,kzL);
				dataLLL_I=flag*src_data[ 2*kxL+1 + nx*kyL+ nxny*kzL ];//get_value_at(2*kxL+1,kyL,kzL);
			    }
			    //      LLU 2
			    kxL= kxLower; kyL=kyLower; int kzU=kzUpper; 
			    float dataLLU_R= 0; float dataLLU_I=0; flag=1;
			    if ( (abs(kxL)<Mid) && (abs(kyL)<Mid) && (abs(kzU)<Mid) ) {//   Step 2 Make sure to be in First BZ
				kxL = (N+kxL)%N;  kyL = (N+kyL)%N; kzU = (N+kzU)%N;
				if (kxL> N/2){kxL=(N-kxL)%N; kyL=(N-kyL)%N ; kzU=(N-kzU)%N ;flag=-1;} // Step 3: use Friedel paired
				dataLLU_R=      src_data[ 2*kxL   + nx*kyL+ nxny*kzU ]; //   image -> get_value_at(2*kxL  ,kyL,kzU);
				dataLLU_I= flag*src_data[ 2*kxL+1 + nx*kyL+ nxny*kzU ]; //   image -> get_value_at(2*kxL+1,kyL,kzU);
			    }
			    //      LUL 3
			    kxL= kxLower; int kyU=kyUpper; kzL=kzLower;  
			    float dataLUL_R= 0; float dataLUL_I=0; flag=1;
			    if ( (abs(kxL)<Mid) && (abs(kyU)<Mid)&& (abs(kzL)<Mid) ) {//  Step 2 Make sure to be in First BZ
				kxL = (N+kxL)%N;  kyU = (N+kyU)%N; kzL = (N+kzL)%N;
				if (kxL> N/2){ kxL=(N-kxL)%N; kyU=(N-kyU)%N; kzL=(N-kzL)%N ;flag=-1;}// Step 3
				dataLUL_R=     src_data[ 2*kxL   + nx*kyU+ nxny*kzL ]; // image -> get_value_at(2*kxL  ,kyU,kzL);
				dataLUL_I=flag*src_data[ 2*kxL+1 + nx*kyU+ nxny*kzL ]; // image -> get_value_at(2*kxL+1,kyU,kzL);
			    }
			    //      LUU 4
			    kxL= kxLower; kyU=kyUpper; kzL=kzUpper;  
			    float dataLUU_R= 0; float dataLUU_I=0; flag=1;
			    if ( (abs(kxL)<Mid) && (abs(kyU)<Mid)&& (abs(kzU)<Mid)) {//   Step 2 Make sure to be in First BZ
				kxL = (N+kxL)%N;  kyU = (N+kyU)%N; kzU = (N+kzU)%N;
				if (kxL> N/2){kxL=(N-kxL)%N; kyU=(N-kyU)%N; kzL=(N-kzL)%N ;flag=-1;} // Step 3
				dataLUU_R=     src_data[ 2*kxL   + nx*kyU+ nxny*kzU ]; // image -> get_value_at(2*kxL  ,kyU,kzU);
				dataLUU_I=flag*src_data[ 2*kxL+1 + nx*kyU+ nxny*kzU ]; // image -> get_value_at(2*kxL+1,kyU,kzU);
			    }
			     //     ULL  5
			    int kxU= kxUpper; kyL=kyLower; kzL=kzLower; 
			    float dataULL_R= 0; float dataULL_I=0; flag=1;
			    if ( (abs(kxU)<Mid) && (abs(kyL)<Mid) && (abs(kzL)<Mid) ) {//    Step 2
				kxU = (N+kxU)%N; kyL = (N+kyL)%N; kzL = (N+kzL)%N;
				if (kxU> N/2){kxU=(N-kxU)%N; kyL=(N-kyL)%N; kzL=(N-kzL)%N ;flag=-1;} // Step 3
				dataULL_R=     src_data[ 2*kxU   + nx*kyL+ nxny*kzL ]; // image -> get_value_at(2*kxU  ,kyL,kzL);
				dataULL_I=flag*src_data[ 2*kxU+1 + nx*kyL+ nxny*kzL ]; // image -> get_value_at(2*kxU+1,kyL,kzL);
			    }
			    //   ULU 6
			    kxU= kxUpper; kyL=kyLower; kzU=kzUpper; 
			    float dataULU_R= 0; float dataULU_I=0; flag=1;
			    if ( (abs(kxU)<Mid) && (abs(kyL)<Mid)&& (abs(kzU)<Mid) ) {//      Step 2
				kxU = (N+kxU)%N; kyL = (N+kyL)%N; kzU = (N+kzU)%N;
				if (kxU> N/2){kxU=(N-kxU)%N; kyL=(N-kyL)%N; kzU=(N-kzU)%N ;flag=-1;} // Step 3
				dataULU_R=     src_data[ 2*kxU   + nx*kyL+ nxny*kzU ]; // image -> get_value_at(2*kxU  ,kyL,kzU);
				dataULU_I=flag*src_data[ 2*kxU+1 + nx*kyL+ nxny*kzU ]; // image -> get_value_at(2*kxU+1,kyL,kzU);
			    }
			    //     UUL 7
			    kxU= kxUpper; kyU=kyUpper; kzL=kzLower;
			    float dataUUL_R= 0; float dataUUL_I=0; flag=1;
			    if ( (abs(kxU)<Mid) && (abs(kyU)<Mid) && (abs(kzL)<Mid) ) {//      Step 2
				kxU = (N+kxU)%N; kyU = (N+kyU)%N; kzL = (N+kzL)%N;
				if (kxU> N/2){kxU=(N-kxU)%N; kyU=(N-kyU)%N; kzL=(N-kzL)%N ;flag=-1;} // Step 3
				dataUUL_R=     src_data[ 2*kxU   + nx*kyU+ nxny*kzL ]; // image -> get_value_at(2*kxU  ,kyU,kzL);
				dataUUL_I=flag*src_data[ 2*kxU+1 + nx*kyU+ nxny*kzL ]; // image -> get_value_at(2*kxU+1,kyU,kzL);
			    }
			    //    UUU 8
			    kxU= kxUpper; kyU=kyUpper; kzU=kzUpper;
			    float dataUUU_R= 0; float dataUUU_I=0; flag=1;
			    if ( (abs(kxU)<Mid) && (abs(kyU)<Mid) && (abs(kzU)<Mid) ) { //       Step 2
				kxU = (N+kxU)%N; kyU = (N+kyU)%N; kzU = (N+kzU)%N;
				if (kxU> N/2) {kxU=(N-kxU)%N; kyU=(N-kyU)%N; kzU=(N-kzU)%N ;flag=-1;} // Step 3
				dataUUU_R=     src_data[ 2*kxU   + nx*kyU+ nxny*kzU ]; // image -> get_value_at(2*kxU  ,kyU,kzU);
				dataUUU_I=flag*src_data[ 2*kxU+1 + nx*kyU+ nxny*kzU ]; // image -> get_value_at(2*kxU+1,kyU,kzU);
			    }
			    //          Step 4    Assign Weights
			    float WLLL = dkxLower*dkyLower*dkzLower ;
			    // WLLL = sqrt(pow(dkxUpper,2)+ pow(dkyUpper,2) + pow(dkzUpper,2)) ;
			    // WLLL = sin(pi* WLLL+.00000001)/(pi* WLLL+.00000001);
			    float WLLU = dkxLower*dkyLower*dkzUpper ;
			    float WLUL = dkxLower*dkyUpper*dkzLower ;
			    float WLUU = dkxLower*dkyUpper*dkzUpper ;
			    float WULL = dkxUpper*dkyLower*dkzLower ;
			    float WULU = dkxUpper*dkyLower*dkzUpper ;
			    float WUUL = dkxUpper*dkyUpper*dkzLower;
			    float WUUU = dkxUpper*dkyUpper*dkzUpper;
			    tempW = WLLL  +  WLLU  +   WLUL  + WLUU  + WULL + WULU +   WUUL +   WUUU ;
			    //           Step 5    Assign Real, then Imaginary Values
			    tempR =  WLLL*dataLLL_R + WLLU*dataLLU_R + WLUL*dataLUL_R + WLUU*dataLUU_R ;
			    tempR += WULL*dataULL_R + WULU*dataULU_R + WUUL*dataUUL_R + WUUU*dataUUU_R ;
			    //
			    tempI  = WLLL*dataLLL_I + WLLU*dataLLU_I + WLUL*dataLUL_I + WLUU*dataLUU_I ;
			    tempI += WULL*dataULL_I + WULU*dataULU_I + WUUL*dataUUL_I + WUUU*dataUUU_I ;
			    //
//			    float phase = -2*pi*(kxNew*xshift+kyNew*yshift+kzNew*zshift)/ny;
			    float tempRb=tempR*Cphase - tempI*Sphase;
			    float tempIb=tempR*Sphase + tempI*Cphase;
			    //
			    des_data[2*kxN    + nx* kyN +nxny*kzN] = tempRb/tempW;
			    des_data[2*kxN+1  + nx* kyN +nxny*kzN] = tempIb/tempW;
		}}}  // end z, y, x loops through new coordinates
	}   //  end  rotations in Fourier Space  3D
	if ((nz > 1)&&(image -> is_complex())&&(zerocorners==5)) {
		//printf("Hello 3-d complex  TransformProcessor \n");
		float phi =  t.get_rotation("eman").get("phi"); phi=pi*phi/180;
		float alt =  t.get_rotation("eman").get("alt"); alt=pi*alt/180;
		float az  =  t.get_rotation("eman").get("az");   az=pi*az /180;
		Vec3f transNow= t.get_trans();
		float xshift= transNow[0]; float yshift= transNow[1];float zshift= transNow[2];
		
		float MatXX = (cos(az)*cos(phi) - sin(az)*cos(alt)*sin(phi) );
		float MatXY = (- cos(az)*sin(phi) - sin(az)*cos(alt)*cos(phi) ) ;
		float MatXZ = sin(az)*sin(alt) ;
		float MatYX = (sin(az)*cos(phi) + cos(az)*cos(alt)*sin(phi) );
		float MatYY = (- sin(az)*sin(phi) + cos(az)*cos(alt)*cos(phi) )  ;
		float MatYZ = - cos(az)*sin(alt) ;
		float MatZX = sin(alt)*sin(phi);
		float MatZY = sin(alt)*cos(phi);
		float MatZZ = cos(alt)  ;
		float tempR=0; float tempI=0; 
		float Mid =(N+1.0)/2.0;  // Check this
		float phaseConstx  = -2*pi*xshift/ny ;
		float k1= cos(phaseConstx); float k2= sin(phaseConstx);
		float k3= 1.0/k1; float k4= k2/k1; // that is 1/cos and tan()
//		float dataRLLL, dataRLLU, dataRLUL, dataRLUU;
//		float dataRULL, dataRULU, dataRUUL, dataRUUU;
//		float dataILLL, dataILLU, dataILUL, dataILUU;
//		float dataIULL, dataIULU, dataIUUL, dataIUUU;
		int nxny = nx*ny;
		
		for (int kzN = 0; kzN < ny; kzN++) { 
		    int kzNew=kzN; 
		    if (kzN >= nx/2) kzNew=kzN-N; //            Step 0 Unalias new coords; 
		    for (int kyN = 0; kyN < ny; kyN++) {//       moves them to lesser mag   
			int kyNew=kyN; 
			if (kyN>=nx/2) kyNew=kyN-ny;  //        Step 0   Unalias
			  //  Step 1: Do inverse Rotation to get former values, and alias   Step1
		        float kxOld =  MatXY * kyNew +  MatXZ *kzNew - MatXX;
	                float kyOld =  MatYY * kyNew +  MatYZ*kzNew  - MatYX;
			float kzOld =  MatZY * kyNew +  MatZZ*kzNew  - MatZX;
    		        float phase = -2*pi*kzNew*zshift/ny-2*pi*kyNew*yshift/ny - phaseConstx ;
			float Cphase = cos(phase);
			float Sphase = sin(phase);
			int kx, ky,kz, II;
			
			int IndexOut= -2+ nx* kyN +nxny*kzN;
			float OutBounds2 = (Mid*Mid- (kyOld*kyOld+kzOld*kzOld)) ;
			
			int kxNewMax= nx/2;
			if (OutBounds2< 0) kxNewMax=0;
			else if (OutBounds2<(nx*nx/4)) kxNewMax=sqrt(OutBounds2);
			for (int kxN = kxNewMax; kxN < nx/2 ; kxN++ ) {
       			      des_data[2*kxN    + nx* kyN +nxny*kzN] = 0;
        		      des_data[2*kxN    + nx* kyN +nxny*kzN+1] = 0;
			}
			
			
			for (int kxNew = 0; kxNew < kxNewMax; kxNew++ ) {
/*			      printf(" kxNew = %d, kyNew = %d, kzNew = %d,kxOld = %3.2f, kyOld = %3.2f, kzOld = %3.2f \n", 
				      kxNew,kyNew,kzNew, kxOld, kyOld, kzOld);*/
			  
			     IndexOut +=2;

			     kxOld +=  MatXX ;
			     kyOld +=  MatYX ;
			     kzOld +=  MatZX ;
                             phase += phaseConstx;
			     Cphase = Cphase*k1 -Sphase*k2; //update using trig addition; this is   cos = cos cos  -sin sin
			     Sphase = Sphase*k3+ Cphase*k4;  //   and   sin = sin  (1/ cos) + cos * tan;
			      
			    //
				if ((abs(kxOld)>=Mid) || (abs(kyOld)>=Mid) || (abs(kzOld)>=Mid) ) { // out of bounds
          			      des_data[IndexOut] = 0;
	        		      des_data[IndexOut+1] = 0;
				   continue;}
/*				if (kxOld*kxOld+OutBounds> 2*Mid*Mid){ // out of bounds
          			      des_data[IndexOut] = 0;
	        		      des_data[IndexOut+1] = 0;
				   continue;}  */
				//
			    int kxLower= floor(kxOld);
			    int kyLower= floor(kyOld);
			    int kzLower= floor(kzOld);
			    //
			    float dkxUpper= (kxOld-kxLower);
			    float dkyUpper= (kyOld-kyLower);
			    float dkzUpper= (kzOld-kzLower);
			    
			    //      LLL  1
			    kx= kxLower; ky =kyLower; kz=kzLower;
//			    dataRLLL=0;   dataILLL=0;
//			    if ( (abs(kx)<Mid) && (abs(ky)<Mid) && (abs(kz)<Mid) ) { //   Step 2 Make sure to be in First BZ
			      int flagLLL=1;
			      if (kx<0) kx += N; if (ky<0) ky += N; if (kz<0) kz += N;
			      if (kx> N/2){kx=N-kx;ky=(N-ky)%N;kz=(N-kz)%N; flagLLL=-1;} // Step 3: if nec, use Friedel paired
			      int kLLL =2*kx   + nx*ky+ nxny*kz ; 
//			      dataRLLL =         src_data[     kLLL ]; 
//			      dataILLL = flagLLL*src_data[ 1 + kLLL ]; //
//			 	} 

			    //      LLU  2
			    kx= kxLower; ky =kyLower; kz=kzLower+1;
//			    dataRLLU=0;	    dataILLU=0;
//			    if ( (abs(kx)<Mid) && (abs(ky)<Mid) && (abs(kz)<Mid) ) { //   Step 2 Make sure to be in First BZ
			      int flagLLU=1;
			      if (kx<0) kx += N; if (ky<0) ky += N; if (kz<0) kz += N;
			      if (kx> N/2){kx=N-kx;ky=(N-ky)%N;kz=(N-kz)%N; flagLLU=-1;} // Step 3: if nec, use Friedel paired
			      int kLLU =2*kx   + nx*ky+ nxny*kz ; 
//			      dataRLLU =         src_data[  kLLU]; 
//			      dataILLU = flagLLU*src_data[1+kLLU];
//			       } 

			    //      LUL  3
			    kx= kxLower; ky =kyLower+1; kz=kzLower;
//			    dataRLUL=0;    dataILUL=0;
//			    if ( (abs(kx)<Mid) && (abs(ky)<Mid) && (abs(kz)<Mid) ) { //   Step 2 Make sure to be in First BZ
			      int flagLUL=1;
			      if (kx<0) kx += N; if (ky<0) ky += N; if (kz<0) kz += N;
			      if (kx> N/2){kx=N-kx;ky=(N-ky)%N;kz=(N-kz)%N; flagLUL=-1;} // Step 3: if nec, use Friedel paired
			      int kLUL =2*kx   + nx*ky+ nxny*kz ; 
//			      dataRLUL =         src_data[  kLUL]; 
//			      dataILUL = flagLUL*src_data[1+kLUL];
//			      } 

			    //      LUU  4
			    kx= kxLower; ky =kyLower+1; kz=kzLower+1;
//			    dataRLUU=0;	    dataILUU=0;
//			    if ( (abs(kx)<Mid) && (abs(ky)<Mid) && (abs(kz)<Mid) ) { //   Step 2 Make sure to be in First BZ
			      int flagLUU=1;
			      if (kx<0) kx += N; if (ky<0) ky += N; if (kz<0) kz += N;
			      if (kx> N/2){kx=N-kx;ky=(N-ky)%N;kz=(N-kz)%N; flagLUU=-1;} // Step 3: if nec, use Friedel paired
			      int kLUU =2*kx   + nx*ky+ nxny*kz ; 
//			      dataRLUU =         src_data[  kLUU]; 
//			      dataILUU = flagLUU*src_data[1+kLUU];
//			      } 

			    //      ULL  5
			    kx= kxLower+1; ky =kyLower; kz=kzLower;
//			    dataRULL=0;	    dataIULL=0;
//			    if ( (abs(kx)<Mid) && (abs(ky)<Mid) && (abs(kz)<Mid) ) { //   Step 2 Make sure to be in First BZ
			      int flagULL=1;
			      if (kx<0) kx += N; if (ky<0) ky += N; if (kz<0) kz += N;
			      if (kx> N/2){kx=N-kx;ky=(N-ky)%N;kz=(N-kz)%N; flagULL=-1;} // Step 3: if nec, use Friedel paired
			      int kULL =2*kx   + nx*ky+ nxny*kz ; 
//			      dataRULL =         src_data[  kULL]; 
//			      dataIULL = flagULL*src_data[1+kULL];
//			      } 

			    //      ULU  6
			    kx= kxLower+1; ky =kyLower; kz=kzLower+1;
//			    dataRULU=0;	    dataIULU=0;
//			    if ( (abs(kx)<Mid) && (abs(ky)<Mid) && (abs(kz)<Mid) ) { //   Step 2 Make sure to be in First BZ
			      int flagULU=1;
			      if (kx<0) kx += N; if (ky<0) ky += N;if (kz<0) kz += N;
			      if (kx> N/2){kx=N-kx;ky=(N-ky)%N;kz=(N-kz)%N; flagULU=-1;} // Step 3: if nec, use Friedel paired
			      int kULU =2*kx   + nx*ky+ nxny*kz ; 
//			      dataRULU =         src_data[  kULU]; 
//			      dataIULU = flagULU*src_data[1+kULU];
//			      } 
		      
			    //      UUL  7
			    kx= kxLower+1; ky =kyLower+1; kz=kzLower;
//			    dataRUUL=0;    dataIUUL=0;
//			    if ( (abs(kx)<Mid) && (abs(ky)<Mid) && (abs(kz)<Mid) ) { //   Step 2 Make sure to be in First BZ
			      int flagUUL=1;
			      if (kx<0) kx += N; if (ky<0) ky += N;if (kz<0) kz += N;
			      if (kx> N/2){kx=N-kx;ky=(N-ky)%N;kz=(N-kz)%N; flagUUL=-1;} // Step 3: if nec, use Friedel paired
			      int kUUL =2*kx   + nx*ky+ nxny*kz ; 
//			      dataRUUL =         src_data[  kUUL]; 
//			      dataIUUL = flagUUL*src_data[1+kUUL];
//			    }

			    //      UUU  8
			    kx= kxLower+1; ky =kyLower+1; kz=kzLower+1;
//			    dataRUUU=0;    dataIUUU=0;
//			    if ( (abs(kx)<Mid) && (abs(ky)<Mid) && (abs(kz)<Mid) ) { //   Step 2 Make sure to be in First BZ
			      int flagUUU=1;
			      if (kx<0) kx += N; if (ky<0) ky += N;if (kz<0) kz += N;
			      if (kx> N/2){kx=N-kx;ky=(N-ky)%N;kz=(N-kz)%N; flagUUU=-1;} // Step 3: if nec, use Friedel paired
			      int kUUU =2*kx   + nx*ky+ nxny*kz ; 
//			      dataRUUU =         src_data[  kUUU]; 
//			      dataIUUU = flagUUU*src_data[1+kUUU];
//			    }

			    //           Step 4    Assign Real, then Imaginary Values
/*			      printf(" kxNew = %d, kyNew = %d, kzNew = %d,kxOld = %3.2f, kyOld = %3.2f, kzOld = %3.2f \n", 
				      kxNew,kyNew,kzNew, kxOld, kyOld, kzOld);
			      printf("  dataRLLL= %f, dataRULL = %f, dataRLUL = %f,dataRUUL = %f, dataRLLU = %f, dataRULU = %f, dataRLUU = %f, dataRUUU = %f \n", 
				      dataRLLL,dataRULL,dataRLUL, dataRUUL,dataRLLU, dataRULU,dataRLUU, dataRUUU);
			      printf("  dataILLL= %f, dataIULL = %f, dataILUL = %f,dataIUUL = %f, dataILLU = %f, dataIULU = %f, dataILUU = %f, dataIUUU = %f \n", 
				      dataILLL,dataIULL,dataILUL, dataIUUL,dataILLU, dataIULU,dataILUU, dataIUUU);*/
/*			      printf("  kLLL= %f, kULL = %d, kLUL = %d,kUUL = %d, kLLU = %d, kULU = %d, kLUU = %d, kUUU = %d \n", 
				      kLLL,kULL,kLUL, kUUL,kLLU, kULU,kLUU, kUUU);*/
			      
			      tempR = Util::trilinear_interpolate(
				 src_data[kLLL] , src_data[kULL], 
				 src_data[kLUL] , src_data[kUUL],
				 src_data[kLLU] , src_data[kULU], 
				 src_data[kLUU] , src_data[kUUU], 
				 dkxUpper,dkyUpper,dkzUpper);
			               
			      tempI = Util::trilinear_interpolate(
				 flagLLL*src_data[kLLL+1], flagULL*src_data[kULL+1], 
				 flagLUL*src_data[kLUL+1], flagUUL*src_data[kUUL+1],
				 flagLLU*src_data[kLLU+1], flagULU*src_data[kULU+1], 
				 flagLUU*src_data[kLUU+1], flagUUU*src_data[kUUU+1], 
				 dkxUpper,dkyUpper,dkzUpper);

			      //           Step 5    Apply translation
			      float tempRb=tempR*Cphase - tempI*Sphase;
			      float tempIb=tempR*Sphase + tempI*Cphase;
			      //
			      des_data[IndexOut]   = tempRb;
			      des_data[IndexOut+1] = tempIb;
		}}}  // end z, y, x loops through new coordinates
	}   //  end  rotations in Fourier Space  3D
	if ((nz > 1)&&(image -> is_real())) {
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
	if(EMData::usecuda == 1 && image->isrodataongpu()){
	        //cout << "using CUDA xform" << endl;
		p = new EMData(0,0,image->get_xsize(),image->get_ysize(),image->get_zsize(),image->get_attr_dict());
		float * m = new float[12];
		Transform inv = t->inverse();
		inv.copy_matrix_into_array(m);
		image->bindcudaarrayA(true);
		p->runcuda(emdata_transform_cuda(m,image->get_xsize(),image->get_ysize(),image->get_zsize()));
		image->unbindcudaarryA();
		delete [] m;
		p->update();
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
	if(EMData::usecuda == 1 && image->isrodataongpu()){
		//cout << "CUDA xform inplace" << endl;
		image->bindcudaarrayA(false);
		float * m = new float[12];
		Transform inv = t->inverse();
		inv.copy_matrix_into_array(m);
		image->runcuda(emdata_transform_cuda(m,image->get_xsize(),image->get_ysize(),image->get_zsize()));
		image->unbindcudaarryA();
		delete [] m;
		use_cpu = false;
		image->update();
	}
#endif
	if ( use_cpu ) {
		float* des_data = transform(image,*t);
		image->set_data(des_data,image->get_xsize(),image->get_ysize(),image->get_zsize());
		image->update();
	}
	float scale = t->get_scale();
	if (scale != 1.0f) {
		image->scale_pixel(1.0f/scale);
//		update_emdata_attributes(image,image->get_attr_dict(),scale);
	}

	if(t) {delete t; t=0;}

	EXITFUNC;
}


void IntTranslateProcessor::assert_valid_aspect(const vector<int>& translation, const EMData* const) const {
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

struct WSsortlist { float pix; short x,y,z; };

int WScmp(const void *a, const void *b) { float x=((WSsortlist*)b)->pix-((WSsortlist*)a)->pix; return (x>0)-(x<0); }		// comparison function for qsort

// inefficient since a copy was probably already made, but best we can do
void WatershedProcessor::process_inplace(EMData *image) {
	EMData *tmp=process(image);
	memcpy(image->get_data(),tmp->get_data(),image->get_size()*sizeof(float));
	delete tmp;
}

EMData *WatershedProcessor::process(const EMData* const image) {
	unsigned int nseg = params.set_default("nseg",12);
	float thr = params.set_default("thr",0.5f);
	int segbymerge = params.set_default("segbymerge",0);
	int verbose = params.set_default("verbose",0);
	if (nseg<=1) throw InvalidValueException(nseg,"nseg must be greater than 1");

	if (segbymerge) { segbymerge=nseg; nseg=4096; }		// set max number of segments to a large (but not infinite) value. Too many will make connectivity matrix too big

	EMData *ret=new EMData(image->get_xsize(),image->get_ysize(),image->get_zsize());
	ret->to_zero();


	int nx=image->get_xsize();
	int ny=image->get_ysize();
	int nz=image->get_zsize();
	if (nz==1) throw ImageDimensionException("Only 3-D data supported");

	// Count the number of above threshold pixels
	size_t n2seg = 0;
	for (int z=0; z<nz; z++) {
		for (int y=0; y<ny; y++) {
			for (int x=0; x<nx; x++) {
				if (image->get_value_at(x,y,z)>=thr) n2seg++;
			}
		}
	}
	if (verbose) printf("%ld voxels above threshold\n",n2seg);

	// Extract the pixels for sorting
//	WSsortlist srt[n2seg];
	WSsortlist * srt = new WSsortlist[n2seg+1];
	size_t i=0;
	for (int z=1; z<nz-1; z++) {
		for (int y=1; y<ny-1; y++) {
			for (int x=1; x<nx-1; x++) {
				if (image->get_value_at(x,y,z)>=thr) {
					srt[i].pix=image->get_value_at(x,y,z);
					srt[i].x=x;
					srt[i].y=y;
					srt[i].z=z;
					i++;
				}
			}
		}
	}
	if (verbose) printf("Voxels extracted, sorting\n");

	// actual sort
	qsort(&srt,n2seg,sizeof(WSsortlist),WScmp);
	if (verbose) printf("Voxels sorted (%1.4g max), starting watershed\n",srt[0].pix);

	// now we start with the highest value and fill in the segments
	float cseg=1.0;
	int start=n2seg;
	for (i=0; i<n2seg; i++) {
		int x=srt[i].x;
		int y=srt[i].y;
		int z=srt[i].z;
		float lvl=0;
		for (int zz=z-1; zz<=z+1; zz++) {
			for (int yy=y-1; yy<=y+1; yy++) {
				for (int xx=x-1; xx<=x+1; xx++) {
					float v=ret->get_value_at(xx,yy,zz);
					if (v>lvl) lvl=v;				// use the highest numbered border segment (arbitrary)
				}
			}
		}
		if (lvl==0) {
			if (verbose) printf("%d %d %d\t%1.0f\t%1.3g\n",x,y,z,cseg,srt[i].pix);
			lvl=cseg;
			cseg+=1.0;
		}
		if (lvl>nseg) {
			start=i;
			if (verbose) printf("Requested number of segments achieved at density %1.4g\n",srt[i].pix);
			break;
		}		// This means we've made as many segments as we need, so we switch to flood-filling
		ret->set_value_at_fast(x,y,z,lvl);
	}

	// We have as many segments as we'll get, but not all voxels have been segmented, so we do a progressive flood fill in density order
	size_t chg=1;
	while (chg) {
		chg=0;
		for (i=start; i<n2seg; i++) {
			int x=srt[i].x;
			int y=srt[i].y;
			int z=srt[i].z;
			if (ret->get_value_at(x,y,z)!=0) continue;	// This voxel is already done

			float lvl=0;
			for (int zz=z-1; zz<=z+1; zz++) {
				for (int yy=y-1; yy<=y+1; yy++) {
					for (int xx=x-1; xx<=x+1; xx++) {
					float v=ret->get_value_at(xx,yy,zz);
					if (v>lvl) lvl=v;				// use the highest numbered border segment (arbitrary)
					}
				}
			}
			if (lvl==0) continue;					// we just skip voxels without any segmented neighbors
			ret->set_value_at_fast(x,y,z,lvl);
			chg+=1;
		}
		if (verbose) printf("%ld voxels changed\n",chg);
	}

	if (segbymerge) {
		if (cseg<segbymerge) return ret;
	}
	else if (cseg<=nseg) return ret;		// We don't have too many segments, so we just return now

	if (verbose) printf("Merging segments\n");
	// If requested, we now merge segments with the most surface contact until we have the correct final number
	if (segbymerge) {
		int nsegstart=(int)cseg;	// number of segments we actually generated
		nseg=(int)cseg;
		EMData *mx=new EMData(nsegstart,nsegstart,1);		// This will be a "contact matrix" among segments
		float *mxd=mx->get_data();

		// each cycle of the while loop, we eliminate one segment by merging
		int sub1=-1,sub2=-1;		// sub2 will be merged into sub1
		nseg++;						// since we don't actually remove one on the first pass, but decrement the counter
		while (segbymerge<nseg) {
			mx->to_zero();

			for (i=0; i<n2seg; i++) {
				int x=srt[i].x;
				int y=srt[i].y;
				int z=srt[i].z;

				int v1=(int)ret->get_value_at(x,y,z);
				if (v1==sub2) { ret->set_value_at_fast(x,y,z,sub1); v1=sub1; }
				mxd[v1+v1*nsegstart]++;					// the diagonal is a count of the number of voxels in the segment
				for (int zz=z-1; zz<=z+1; zz++) {
					for (int yy=y-1; yy<=y+1; yy++) {
						for (int xx=x-1; xx<=x+1; xx++) {
							int v2=(int)ret->get_value_at(xx,yy,zz);
							if (v2==sub2) v2=sub1;		// pretend that any sub2 values are actually sub1
							if (v1==v2) continue;
							mxd[v1+v2*nsegstart]+=image->get_value_at(xx,yy,zz);		// We weight the connectivity by the image value
						}
					}
				}
			}
			mx->update();
			nseg--;					// number of segments left
			if (verbose && sub1==-1) { mx->write_image("contactmx.hdf",0); } 		// for debugging

			sub1=-1;
			sub2=-1;
			// contact matrix complete, now figure out which 2 segments to merge
			// diagonal of matrix is a count of the 'volume' of the segment. off-diagonal elements are surface area of contact region (roughly)
			// we want to normalize the surface area elements so they are roughly proportional to the size of the segment, so we don't merge
			// based on total contact area, but contact area as a fraction of the total area.
			float bestv=-1.0;
			for (int s1=1; s1<nsegstart; s1++) {
				for (int s2=1; s2<nsegstart; s2++) {
					if (s1==s2) continue;				// ignore the diagonal
					float v=mxd[s1+s2*nsegstart];
					if (v==0) continue;					// empty segment
//					v/=(pow(mxd[s1+s1*nsegstart],0.6667f)+pow(mxd[s2+s2*nsegstart],0.6667f));	// normalize by the sum of the estimated surface areas (no shape effects)
					v/=max(mxd[s1+s1*nsegstart],mxd[s2+s2*nsegstart]);	// normalize by the sum of the estimated surface areas (no shape effects)
					if (v>bestv) { bestv=v; sub1=s1; sub2=s2; }
			 	}
			}
			float mv=0;
			int mvl=0;
			for (i=nsegstart+1; i<nsegstart*nsegstart; i+=nsegstart+1) 
				if (mxd[i]>mv) { mv=mxd[i]; mvl=i/nsegstart; }
			if (verbose) printf("Merging %d to %d (%1.0f, %d)\n",sub2,sub1,mv,mvl);
			if (sub1==-1) {
				if (verbose) printf("Unable to find segments to merge, aborting\n");
				break;
			}
		}

	}

	delete [] srt;

	return ret;
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
	if (EMData::usecuda == 1 && image->getcudarwdata()) {
		//cout << "CUDA rotate 180" << endl;
		emdata_rotate_180(image->getcudarwdata(), image->get_xsize(), image->get_ysize());
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
	bool tozero = params.set_default("tozero",false);
	float new_min_vals = min;
	float new_max_vals = max;
	if (tomean) new_min_vals = new_max_vals = image->get_attr("mean");
	if (tozero) new_min_vals = new_max_vals = 0.0;

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
	params.put("minval",mean - nsigma*sigma);
	params.put("maxval",mean + nsigma*sigma);

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
	image->process_inplace("xform.phaseorigin.tocenter");

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
	image->process_inplace("xform.phaseorigin.tocenter");

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
	float len = params.set_default("length", 16.2f); //in angstroms
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
	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i, ++dat) {
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

	for (int k = 0; k < nz; ++k) //taking slices along z axis
	{
		rho_x_sum = rho_y_sum = rho_sum = 0; //Set to zero for a new slice

		if (k >= z_start && k <= z_stop)
		//Apply the radial profile only between z_start and z_stop on the z axis
		{
			//Calculating CM for the slice...
			for (int j = 0; j < ny; ++j)
			{
				for (int i = 0; i < nx; ++i, ++dat)
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
				for (int j=0; j<ny;++j)
				{
					for (int i=0;i<nx;++i,++dat)
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

//Do NOT use this code as an example, (const EMData* const image) needs to be used otherwise process_inplace is called and the params are throw out!
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

EMData* ConvolutionKernelProcessor::process(const EMData* const image)
{
	if (image->get_zsize()!=1) throw ImageDimensionException("Only 2-D images supported");

	EMData* conv = new EMData(image->get_xsize(),image->get_ysize(),1);
	vector<float>kernel = params["kernel"];

	if (fmod(sqrt((float)kernel.size()), 1.0f) != 0) throw InvalidParameterException("Convolution kernel must be square!!");

	float* data = image->get_data();
	float* cdata = conv->get_data();	// Yes I could use set_value_at_fast, but is still slower than this....

	//I could do the edges by wrapping around, but this is not necessary(such functionality can be iplemented later)
	int ks = int(sqrt(float(kernel.size())));
	int n = (ks - 1)/2;
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	for (int i = n; i < (nx - n); i++) {
		for (int j = n; j < (ny - n); j++) {
			//now do the convolution
			float cpixel = 0;
			int idx = 0;
			// Perahps I could use some ofrm of Caching to speed things up?
			for (int cx = -n; cx <= n; cx++) {
				for (int cy = -n; cy <= n; cy++) {
					cpixel += data[(i+cx) + (j+cy) * nx]*kernel[idx];
					idx++;
				}
			}
			cdata[i + j * nx] = cpixel;
		}
	}

	return conv;
}

void ConvolutionKernelProcessor::process_inplace(EMData * image )
{
	throw UnexpectedBehaviorException("Not implemented yet");

	return;
}


EMData* RotateInFSProcessor::process(const EMData* const image) //
{

    EMData* imageCp        = image -> copy(); // This is the rotated image
    process_inplace(imageCp);

    return imageCp;
}


void RotateInFSProcessor::process_inplace(EMData * image) // right now for 2d images
{
//	float angle = params["angle"];


//	Transform *rotNow  = params.set_default("transform",&Transform());
	Transform *rotNow  = params["transform"];
	float interpCutoff = params.set_default("interpCutoff",0.8f);   // JFF, can you move this to input parameter?
//	float interpCutoff = params["interpCutoff"];   // JFF, can you move this to input parameter?
//	float interpCutoff =.8;   // JFF, can you move this to input parameter?


	int debug=0;

//	error: conversion from 'EMAN::EMObject' to non-scalar type 'EMAN::Transform' requested


	// if 2N is size of image, then sizes of FFT are (2N+2,2N,2N) or  (2N+2,2N,1)
	// if 2N+1 is size of image, then sizes of FFT are (2N+2,2N+1,2N+1) or  (2N+2,2N+1,1)

	int x_size = image->get_xsize();  //16
	int y_size = image->get_ysize();  int y_odd=  (y_size%2);
	int z_size = image->get_zsize();

//	float size_check = abs(y_size-z_size)+abs(x_size-y_size-2+y_odd);
//	if (size_check!=0) throw ImageDimensionException("Only cubic images");
	int size_check = abs(x_size-y_size-2+y_odd)+ abs(z_size-1)*abs(z_size-y_size);
	int N = x_size/2-1;
	int Mid = N+1;
	if (size_check!=0) throw ImageDimensionException("Only square or cubic  images for now");
	if (image->is_real()) throw ImageDimensionException("Only for Fourier images");
//	if (y_odd==0) throw ImageDimensionException("Only use odd images for now");

	if (debug) printf("Mid=%d, x_size=%d, y_size=%d, N=%d, z_size=%d \n", Mid,x_size,y_size, N, z_size );

	EMData* RotIm        = image -> copy(); // This is the rotated image
	EMData* WeightIm     = image -> copy(); WeightIm     ->to_zero();// This is the Weight image for rotated image
	EMData* PhaseOrigIm  = new EMData(N+1,2*N+1,z_size) ; PhaseOrigIm  ->to_zero();// This is the Weight image for rotated image
	EMData* PhaseFinalIm = PhaseOrigIm -> copy(); PhaseFinalIm ->to_zero();// This is the Weight image for rotated image
	EMData* MagFinalIm   = PhaseOrigIm -> copy(); MagFinalIm   ->to_zero();// This is the Weight image for rotated image

//	float*   data      = image    -> get_data();
//	float* Rotdata     = RotIm    -> get_data();  // This is the data of the rotated image
//	float* WeightData  = WeightIm -> get_data();  //

	float  WeightNowX, WeightNowY, WeightNowZ ;
	int    kxMin,kxMax, kyMin, kyMax,  kzMin, kzMax, kxBefore, kyBefore, kzBefore ;
	float  kxRT, kyRT, kzRT ;

	Vec3f PosAfter;
	Vec3f PosBefore;
	Transform invRotNow;
	// Fill out real and imaginary full images

	//int kz=0;

	if (debug) {image -> write_image("OrigImageFFT.hdf");
	            printf("Just wrote OrigImageFFT.hdf \n"); }


	for (kxBefore = 0; kxBefore <=  N      ; ++kxBefore) {  // These are the  kx, ky coordinates of the original image
	   for (kyBefore = 0; kyBefore < y_size ; ++kyBefore)  {         //  We need to rephase
	      for (kzBefore = 0; kzBefore < z_size ; ++kzBefore)  {         //  We need to rephase

	    // Now we need a
	      float CurrentReal = RotIm -> get_value_at(2*kxBefore  ,kyBefore, kzBefore);
	      float CurrentImag = RotIm -> get_value_at(2*kxBefore+1,kyBefore, kzBefore);

//         fftOfxPRB3(1+mod(ik-Mid,N))=Grand*exp(-2*pi*1i*Mid*(ik-Mid)/x_size); % Phase to apply to centered version

//	      float Phase    = -2*pi*(kxBefore+kyBefore + kzBefore)*(Mid)/y_size;
	      float Phase    = -pi*(kxBefore+kyBefore + kzBefore)*x_size/y_size;
	      // Phase    = 0;
	      float CosPhase = cos( Phase);
	      float SinPhase = sin( Phase);

	      float NewRealValue = CosPhase*CurrentReal -SinPhase*CurrentImag;
	      float NewImagValue = SinPhase*CurrentReal +CosPhase*CurrentImag;

	      RotIm ->set_value_at(2*kxBefore  ,kyBefore, kzBefore, NewRealValue);
	      RotIm ->set_value_at(2*kxBefore+1,kyBefore, kzBefore, NewImagValue);
	}}}

	if (debug) {RotIm  -> write_image("OrigImageFFTAfterPhaseCorrection.hdf");
	            printf("  Just wrote OrigImageFFTAfterPhaseCorrection.hdf \n");}

        // RotIm ->set_value_at(2*Mid-1,0, 0, 0);
        if (debug) printf("  Just about to start second loop  \n");

	image ->to_zero();
        invRotNow = rotNow ->inverse(); //  no match for 'operator=' in 'PosBefore = EMAN::operator*(const EMAN::Transform&, const EMAN::Transform&)((


	for (int kxAfter = 0; kxAfter <= N  ; ++kxAfter) {  // These are the  kx, ky, kz coordinates of the rotated image
	  for (int kyAfter = -N; kyAfter < y_size-N     ; ++kyAfter)  {     // referring to a properly centered version
	     for (int kzAfter = -z_size/2; kzAfter <= z_size/2  ; ++kzAfter)  {

	    // Now we need a

	      PosAfter = Vec3f(kxAfter, kyAfter, kzAfter);
	      PosBefore = invRotNow*PosAfter;
	      kxRT = PosBefore[0]; // This will be the off-lattice site, where the point was rotated from
	      kyRT = PosBefore[1]; //
	      kzRT = PosBefore[2]; //


	      kxMin = ceil( kxRT-interpCutoff); kxMax = floor(kxRT+interpCutoff);
	      kyMin = ceil( kyRT-interpCutoff); kyMax = floor(kyRT+interpCutoff);
	      kzMin = ceil( kzRT-interpCutoff); kzMax = floor(kzRT+interpCutoff);


//              printf("Block 0,kx=%d, ky=%d,kxMin=%d, kyMin=%d, kxMax=%d, kyMax=%d, kyAfter=%d  \n",kxAfter,kyAfter,kxMin,kyMin, kxMax, kyMax, kyAfter);
	      //continue;
	      for (int kxI= kxMin; kxI <= kxMax; ++kxI){  // go through this
		for (int kyI= kyMin; kyI <= kyMax; ++kyI){  // and get values to interp
		   for (int kzI= kzMin; kzI <= kzMax; ++kzI){ //

//
		     if (abs(kxI) >N) continue; // don't go off the lattice
		     if (abs(kyI) >N) continue;
		     if (abs(kzI) >z_size/2) continue;

		     float distx= abs(kxI-kxRT);
		     float disty= abs(kyI-kyRT);
		     float distz= abs(kzI-kzRT);

		     // fold kxI, kyI back into lattice if possible
		     int IsComplexConj= 1;

		     if (kxI<0) IsComplexConj=-1;
		     kxBefore= IsComplexConj*kxI; // The Proper coordinates will be
		     kyBefore= IsComplexConj*kyI; // where the original data is stored
		     kzBefore= IsComplexConj*kzI; // At this point kxBefore >=0, but not necessarily kyBefore

		     if ( kyBefore<0 ) kyBefore += y_size; // makes sure kyBefore is also >0
		     if ( kzBefore<0 ) kzBefore += y_size; // makes sure kzBefore is also >0

		     WeightNowX  = (distx ==0)? 1: (sin(pi*distx) /(pi*distx)) ;
		     WeightNowY  = (disty ==0)? 1: (sin(pi*disty) /(pi*disty)) ;
		     WeightNowZ  = (distz ==0)? 1: (sin(pi*distz) /(pi*distz)) ;


		     float CurrentValue;
		     float ToAdd;
		     int kyAfterInd = (kyAfter+y_size)%(y_size);
		     int kzAfterInd = (kzAfter+z_size)%(z_size);

		    // if (kxAfter==0) IsComplexConj*=-1;

//		     if ((kxI+kyI)%1 ==0)
//		         printf("Block5,kx=%d, ky=%d,kxI=%d, kyI=%d, kxBefore=%d, kyBefore=%d  \n",kxAfter,kyAfter,kxI,kyI, kxBefore,kyBefore);
//		         printf("  %d,     %d,  %d,  %d,        %d,  %d,       %d, %d \n",IsComplexConj,kxAfter,kyAfter, kyAfterInd,kxI,kyI, kxBefore,kyBefore);

		     CurrentValue =   image -> get_value_at(2*kxAfter,kyAfterInd, kzAfterInd);  // Update real part of Image
		     ToAdd =   WeightNowX*WeightNowY*WeightNowZ*(RotIm -> get_value_at(2*kxBefore,kyBefore, kzBefore));
		     image -> set_value_at(2*kxAfter  ,kyAfterInd  , kzAfterInd,  ToAdd   + CurrentValue );


		     CurrentValue =   WeightIm -> get_value_at(kxAfter,kyAfterInd, kzAfterInd);    // Update real  part of Weight image
		     ToAdd =   WeightNowX*WeightNowY;
		     WeightIm -> set_value_at(kxAfter  , kyAfterInd , kzAfterInd,  abs(ToAdd)   + CurrentValue );

		     CurrentValue = image -> get_value_at(2*kxAfter+1,kyAfterInd);    // Update imaginary   part of image
		     ToAdd =   IsComplexConj*WeightNowX*WeightNowY*RotIm -> get_value_at(2*kxBefore+1,kyBefore, kzBefore );
		     image -> set_value_at(2*kxAfter+1  , kyAfterInd , kzAfterInd,  ToAdd   + CurrentValue );


		}}}


	  }}}

//        Set the image values to the rotated image, because we do it in place


	if (debug) { image        -> write_image("RotImageBeforeFinalPhaseCorrection.hdf");
	             printf("  Just wrote RotImageBeforeFinalPhaseCorrection.hdf \n");   }


	for (kxBefore = 0; kxBefore <= N     ; ++kxBefore) {      // This is  the normalization step
	  for (kyBefore = 0; kyBefore < y_size   ; ++kyBefore)  {  // These are the  kx, ky, kz coordinates of the original image
              for (kzBefore = 0; kzBefore < z_size   ; ++kzBefore)  {

	      float CurrentReal = image -> get_value_at(2*kxBefore   , kyBefore, kzBefore);
	      float CurrentImag = image -> get_value_at(2*kxBefore+1 , kyBefore, kzBefore);

              PhaseFinalIm -> set_value_at(kxBefore,kyBefore, kzBefore, atan2(CurrentImag,CurrentReal));
              MagFinalIm   -> set_value_at(kxBefore,kyBefore, kzBefore, sqrt(CurrentImag*CurrentImag+CurrentReal*CurrentReal) );
	      float WeightNow    = WeightIm -> get_value_at(kxBefore,kyBefore, kzBefore);
	      if (WeightNow>0) {
		 float val =  (image->get_value_at(2*kxBefore   , kyBefore, kzBefore))/WeightNow;
	         image -> set_value_at(2*kxBefore   , kyBefore, kzBefore, val);
		 val      =  (image->get_value_at(2*kxBefore +1  , kyBefore, kzBefore))/WeightNow;
	         image -> set_value_at(2*kxBefore +1  , kyBefore, kzBefore,  val);
	      }

	}}}

	if (debug) { printf("  Just did normalization step \n");}


	for ( kxBefore = 0; kxBefore < Mid     ; ++kxBefore) {     //  This is the rephase step
	  for ( kyBefore = 0; kyBefore < y_size   ; ++kyBefore)  {
  	    for ( kzBefore = 0; kzBefore < z_size   ; ++kzBefore)  {

	      float CurrentReal = image -> get_value_at(2*kxBefore  ,kyBefore, kzBefore);
	      float CurrentImag = image -> get_value_at(2*kxBefore+1,kyBefore, kzBefore);

//	      float Phase    = +2*pi*(kxBefore+kyBefore+kzBefore)*(Mid)/y_size;
	      float Phase    =  pi*(kxBefore + kyBefore + kzBefore)*x_size/y_size;
	      // Phase    = 0; // Offset should be Mid-1
	      float CosPhase = cos( Phase);
	      float SinPhase = sin( Phase);

	      float NewRealValue = CosPhase*CurrentReal -SinPhase*CurrentImag;
	      float NewImagValue = SinPhase*CurrentReal +CosPhase*CurrentImag;

	      image ->set_value_at(2*kxBefore,  kyBefore,  kzBefore, NewRealValue);
	      image ->set_value_at(2*kxBefore+1,kyBefore,  kzBefore, NewImagValue);
	}}}

	if (debug) {
	   image        -> write_image("RotatedImageFFT.hdf");
	   PhaseFinalIm -> write_image("PhaseImInFS.hdf");   // These are the phases,mags of the image when properly centered
	   MagFinalIm   -> write_image("MagFinalImInFS.hdf");
	   WeightIm     -> write_image("WeightIm.hdf");
	   printf("  Just wrote RotatedImageFFT.hdf \n");
	}

	image -> update();

}
	
EMData* CircularAverageBinarizeProcessor::process(const EMData* const image)  // now do it layer by layer in 3d
{
	int thr=params.set_default("thresh",5);
	
	EMData* bwmap= image -> copy(); 
	int x_size = image->get_xsize(); 
	int y_size = image->get_ysize();
	int z_size = image->get_zsize();

	int ix,iy,iz,it,count,ic;
	int *dx=new int[thr*8],*dy=new int[thr*8];
	for (it=1; it<=thr; it++){
		// calculate the indexes
		count=0;
		for (ix=-thr-1; ix<=thr+1; ix++){
			for (iy=-thr-1; iy<=thr+1; iy++){
				int d2=ix*ix+iy*iy;
				if (d2>=it*it && d2<(it+1)*(it+1)){
					dx[count]=ix;
					dy[count]=iy;
					count++;
				}
			}
		}
		// for each pixel
		for (iz=0; iz<z_size; iz++){
			for (ix=0; ix<x_size; ix++){
				for (iy=0; iy<y_size; iy++){
					// circular average on each ring
					float mn=0;
					if (bwmap->get_value_at(ix,iy,iz)==0)
						continue;
					for (ic=0; ic<count; ic++){
						mn+=image->sget_value_at(ix+dx[ic],iy+dy[ic],iz);
					}
					mn/=count;
					if (mn>bwmap->get_value_at(ix,iy,iz))
						mn=0;
					bwmap->set_value_at(ix,iy,iz,mn);
				}
			}
		}
	}
	// binarize image
	for (iz=0; iz<z_size; iz++){
		for (ix=0; ix<x_size; ix++){
			for (iy=0; iy<y_size; iy++){
				if (bwmap->get_value_at(ix,iy,iz)>0)
					bwmap->set_value_at(ix,iy,iz,1);
				else
					bwmap->set_value_at(ix,iy,iz,0);
			}
		}
	}
	delete 
	dx;
	delete dy;
	return bwmap;
	
	
}
void CircularAverageBinarizeProcessor::process_inplace(EMData * image)
{
	EMData *tmp=process(image);
	memcpy(image->get_data(),tmp->get_data(),(size_t)image->get_xsize()*image->get_ysize()*image->get_zsize()*sizeof(float));
	delete tmp;
	image->update();
	return;
	
}

EMData* ObjDensityProcessor::process(const EMData* const image) //
{

    EMData* imageCp= image -> copy(); 
    process_inplace(imageCp);

    return imageCp;
}

void ObjDensityProcessor::process_inplace(EMData * image)
{
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	float threshold=params.set_default("thresh",1.5);
	bool nb8=params.set_default("more_neighbor",false);
	
	// label each object first
	EMData *label=image->process("morph.object.label",Dict("thresh",threshold,"more_neighbor",nb8));
	int nobj=int(label->get_attr("maximum"))+1;
	float *sden=new float[nobj];	// sum density of each object
	for (int i=0; i<nobj; i++) sden[i]=0;
	
	for (int x=0; x<nx; x++){
		for (int y=0; y<ny; y++){
			for (int z=0; z<nz; z++){
				float v=image->get_value_at(x,y,z);
				if (v<threshold)
					continue;
				
				int li=label->get_value_at(x,y,z);
				sden[li]+=v;
			}
		}
	}
	
	for (int x=0; x<nx; x++){
		for (int y=0; y<ny; y++){
			for (int z=0; z<nz; z++){
				int li=label->get_value_at(x,y,z);
				if (li==0)
					continue;
				image->set_value_at_fast(x,y,z,sden[li]);
			}
		}
	}
	delete label;
	delete []sden;
}
EMData* ObjLabelProcessor::process(const EMData* const image) //
{

    EMData* imageCp= image -> copy(); 
    process_inplace(imageCp);

    return imageCp;
}

void ObjLabelProcessor::process_inplace(EMData * image)
{
	// This structure is copied from AutoMaskDustProcessor
	
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	float threshold=params.set_default("thresh",1.5);
	bool nb8=params.set_default("more_neighbor",false);
	bool writecenter=params.set_default("write_centers",false);

	vector<float> centers;
	EMData *mask = new EMData();
	mask->set_size(nx, ny, nz);
	mask->to_one();
	int count=0;
	int zs=(nz>1)? 1 : 0;
	for (int zz = 0; zz < nz; zz++) {
		for (int yy = 0; yy < ny; yy++) {
			for (int xx = 0; xx < nx; xx++) {
				if (image->get_value_at(xx,yy,zz)>threshold && mask->get_value_at(xx,yy,zz)==1.0) {
					count++;
					vector<Vec3i> pvec;
					pvec.push_back(Vec3i(xx,yy,zz));
					for (uint i=0; i<pvec.size(); i++) {
						// Duplicates will occur the way the algorithm is constructed, so we eliminate them as we encounter them
						if (mask->sget_value_at(pvec[i])==0.0f) {
							pvec.erase(pvec.begin()+i);
							i--;
							continue;
						}

						// mask out the points in the volume
						mask->set_value_at(pvec[i],0.0f);

						int x=pvec[i][0];
						int y=pvec[i][1];
						int z=pvec[i][2];
						// Any neighboring values above threshold we haven't already set get added to the list
						if (nb8){
							for (int ix=x-1; ix<=x+1; ix++){
								for (int iy=y-1; iy<=y+1; iy++){
									for (int iz=z-zs; iz<=z+zs; iz++){
										if (image->sget_value_at(ix,iy,iz)>threshold && mask->sget_value_at(ix,iy,iz)==1.0){
											pvec.push_back(Vec3i(ix,iy,iz));
										}
									}
									
								}
							}
						}
						else{
							if (image->sget_value_at(x-1,y,z)>threshold && mask->sget_value_at(x-1,y,z)==1.0){
								pvec.push_back(Vec3i(x-1,y,z));
							}
							if (image->sget_value_at(x+1,y,z)>threshold && mask->sget_value_at(x+1,y,z)==1.0){
								pvec.push_back(Vec3i(x+1,y,z));
							}							
							if (image->sget_value_at(x,y-1,z)>threshold && mask->sget_value_at(x,y-1,z)==1.0){
								pvec.push_back(Vec3i(x,y-1,z));
							}
							if (image->sget_value_at(x,y+1,z)>threshold && mask->sget_value_at(x,y+1,z)==1.0){
								pvec.push_back(Vec3i(x,y+1,z));
							}
							if (image->sget_value_at(x,y,z-1)>threshold && mask->sget_value_at(x,y,z-1)==1.0){
								pvec.push_back(Vec3i(x,y,z-1));
							}
							if (image->sget_value_at(x,y,z+1)>threshold && mask->sget_value_at(x,y,z+1)==1.0){ 
								pvec.push_back(Vec3i(x,y,z+1));
							}
						}
					}

					for (uint i=0; i<pvec.size(); i++) mask->set_value_at(pvec[i],2.0);
					for (uint i=0; i<pvec.size(); i++) image->set_value_at(pvec[i],count);
					if (writecenter){
						float cnt[3]={0,0,0};
						for (uint i=0; i<pvec.size(); i++){
	// 						printf("%d,%d,%d\n",pvec[i][0],pvec[i][1],pvec[i][2]);
							cnt[0]+=pvec[i][0];
							cnt[1]+=pvec[i][1];
							cnt[2]+=pvec[i][2];
						}
						cnt[0]/=pvec.size();
						cnt[1]/=pvec.size();
						cnt[2]/=pvec.size();
	// 					printf("%f,%f,%f\n",cnt[0],cnt[1],cnt[2]);
						centers.push_back(cnt[0]);
						centers.push_back(cnt[1]);
						centers.push_back(cnt[2]);
					}
				}
			}
		}
	}
	if (writecenter)
		image->set_attr("obj_centers",centers);

	delete mask;
}

EMData* BwThinningProcessor::process(const EMData* const image) //
{

    EMData* imageCp= image -> copy(); 
    process_inplace(imageCp);

    return imageCp;
}

void BwThinningProcessor::process_inplace(EMData * image){

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	
	float threshold=params.set_default("thresh",0);
	int ntimes=params.set_default("ntimes",-1);
	int verbose=params.set_default("verbose",0);
	int preserve=params.set_default("preserve_value",false);
	if (nz > 1) {
		ImageDimensionException("Only 2-D images supported");
	}
	float array[9];
	int n=1;
	EMData* imageCp;
	if (preserve)
		imageCp= image -> copy(); 
	float *data = image->get_data();
	size_t total_size = (size_t)nx * (size_t)ny * (size_t)nz;
	// binarize image and deal with boundary points first
	for (int j=0; j<ny; j++){
		int jnx = j * nx;
		for (int i=0; i<nx; i++){
			if(i==0 || i==nx-1 || j==0 || j==ny-1)
				data[i+jnx]=0;
			else{
				if (data[i+jnx]>threshold)
					data[i+jnx]=1;
				else
					data[i+jnx]=0;
			}
		}
	}
	
	float *data2 = new float[total_size];
	int ntstep=1,allt=ntimes;
	if (ntimes<0){	// thin to skeleton
		allt=1;
		ntstep=0;
	}
	// thinning
	int cg;
	for (int nt=0; nt<allt; nt+=ntstep){
		cg=0;
		for (int st = 0; st<2; st++){
			memcpy(data2, data, total_size * sizeof(float));
			for (int j = n; j < ny - n; j++) {
				int jnx = j * nx;
				for (int i = n; i < nx - n; i++) {
					size_t s = 0;
					for (int i2 = i - n; i2 <= i + n; i2++) {
						for (int j2 = j - n; j2 <= j + n; j2++) {
							array[s] = data2[i2 + j2 * nx];
							++s;
						}
					}

					cg+=process_pixel(&data[i + jnx ], array, st);
				}
			}
		}
		if (verbose>0)
			printf("%d pixels changed\n",cg);
		if(cg==0)
			break;
	}
	
	// remove corner pixels when doing skeletonization
	if (ntimes<0){
		cg=0;
		memcpy(data2, data, total_size * sizeof(float));
		for (int j = n; j < ny - n; j++) {
			int jnx = j * nx;
			for (int i = n; i < nx - n; i++) {
				size_t s = 0;
				for (int i2 = i - n; i2 <= i + n; i2++) {
					for (int j2 = j - n; j2 <= j + n; j2++) {
						array[s] = data2[i2 + j2 * nx];
						++s;
					}
				}

				cg+=process_pixel(&data[i + jnx ], array, 2);
			}
		}
	}
	if (verbose>0)
		printf("%d corner pixels\n",cg);

	image->update();

	if( data2 )
	{
		delete[]data2;
		data2 = 0;
	}
	
	if (preserve){
		image->mult(*imageCp);
		delete imageCp;
	}
	
}


int BwThinningProcessor::process_pixel(float* data, float* array, int step){
	
	if (*data==0){
		return 0;
	}
	int bp=-1; // number of 1 neighbors, not counting itself 
	for (int i=0; i<9; i++){
		if (array[i]>0)
			bp++;
	}
	if (bp<2 || bp>6){
		return 0;
	}
	int ap=0; // number of transitions from 0 to 1
	int order[9]={0,1,2,5,8,7,6,3,0};
	for (int i=0; i<8; i++){
		if (array[order[i]]==0 && array[order[i+1]]>0){
			ap++;
		}
	}
	if (ap!=1 && step<2)
		return 0;
	
	if (step==0){
		if(array[order[1]]*array[order[3]]*array[order[5]]>0)
			return 0;
		if(array[order[3]]*array[order[5]]*array[order[7]]>0)
			return 0;		
	}
	
	if (step==1){
		if(array[order[1]]*array[order[3]]*array[order[7]]>0)
			return 0;
		if(array[order[1]]*array[order[5]]*array[order[7]]>0)
			return 0;		
	}
	
	if (step==2){
		if (bp==2){
			if(array[order[1]]*array[order[3]]>0 
			|| array[order[3]]*array[order[5]]>0
			|| array[order[5]]*array[order[7]]>0
			|| array[order[7]]*array[order[1]]>0
			){
				*data=0;
				return 1;
			}
			else{
				return 0;
			}
		}
		if (ap==2 ){
			if(array[order[1]]*array[order[3]]>0 
			|| array[order[3]]*array[order[5]]>0
			){
				*data=0;
				return 1;
			}
			else{
				return 0;
			}
		}
		return 0;
	}
	
	*data=0;
	return 1;
	
}

EMData* PruneSkeletonProcessor::process(const EMData* const image) //
{
	EMData* imageCp= image -> copy(); 
	process_inplace(imageCp);
	return imageCp;
}

void PruneSkeletonProcessor::process_inplace(EMData * image){
	// This function is far from optimized, but it works...
	
	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();
	
	float threshold=params.set_default("thresh",0);
	int verbose=params.set_default("verbose",0);
	int maxdist=params.set_default("maxdist",3);
	
	if (nz > 1) {
		ImageDimensionException("Only 2-D images supported");
	}
	
	float *data = image->get_data();
	size_t total_size = (size_t)nx * (size_t)ny * (size_t)nz;
	
	float *data2 = new float[total_size];
	memcpy(data2, data, total_size * sizeof(float));
	
	// binarize image first
	for (int j=0; j<ny; j++){
		int jnx = j * nx;
		for (int i=0; i<nx; i++){
		
			if (data[i+jnx]>threshold)
				data[i+jnx]=1;
			else
				data[i+jnx]=0;
		}
	}
	
	float array[9];
	image->to_zero();
	
	// find branch points
	int nbranch=0;
	for (int j=1; j < ny-1; j++) {
		int jnx = j * nx;
		for (int i=1; i<nx-1; i++) {
			
			if (data2[i+jnx]<=threshold)
				continue;
			int s=0;
			for (int i2 = i-1; i2 <= i + 1; i2++) {
				for (int j2 = j - 1; j2 <= j + 1; j2++) {
					array[s++] = data2[i2 + j2 * nx];
				}
			}
			int ap=0; // number of transitions from 0 to 1
			int order[9]={0,1,2,5,8,7,6,3,0};
			for (int oi=0; oi<8; oi++){
				if (array[order[oi]]<=threshold && array[order[oi+1]]>threshold){
					ap++;
				}
			}
			if (ap>2){
				data[i+jnx]=1;
				nbranch++;
			}
		}
	}
	if (verbose>0)
		printf("%d branch pixels\n",nbranch);
	
	// now, data->branch points, data2->binarized image
	
	// distance to branch point
	for (int j=1; j<ny; j++){
		int jnx=j*nx;
		for (int i=0; i<nx; i++){
			data[i+jnx]=(1-data[i+jnx])*(maxdist+1);
		}
		
	}
	for (int dt=0; dt<maxdist; dt++){
		for (int j=1; j < ny-1; j++) {
			int jnx=j*nx;
			for (int i=1; i<nx-1; i++) {
				
				if (data2[i+jnx]<=threshold)
					continue;
				if (data[i+jnx]<=maxdist)
					continue;
				int db=maxdist;	// distance from nearest branch point. 
				for (int i2=i-1; i2<=i+1; i2++) {
					for (int j2=j-1; j2<=j+1; j2++) {
						db=(data[i2+j2*nx]==dt ? dt : db);
					}
				}
				data[i+jnx]=db+1;
			}
		}
	}
	// now, data->distance to the nearest branch point
	
	// mark endpoints for deletion
	int nend=0;
	for (int j=1; j < ny-1; j++) {
		int jnx=j*nx;
		for (int i=1; i<nx-1; i++) {
			
			if (data2[i+jnx]<=threshold)
				continue;
			if (data[i+jnx]>maxdist)
				continue;
			int nb=-1;	// number of neighbors
			for (int i2=i-1; i2<=i+1; i2++) {
				for (int j2=j-1; j2<=j+1; j2++) {
					nb+=(data2[i2+j2*nx]>threshold ? 1 : 0);
				}
			}
			if (nb==1){	// endpoint found
				data[i+jnx]=-data[i+jnx];	//mark for deletion
				data2[i+jnx]=threshold;
				nend++;
			}
		}
	}
	
	// remove marked branches
	for (int dt=-maxdist; dt<-1; dt++){
		for (int j=1; j < ny-1; j++) {
			int jnx=j*nx;
			for (int i=1; i<nx-1; i++) {
				
				if (data2[i+jnx]<=threshold)
					continue;
				if (data[i+jnx]<=0)
					continue;
				int rm=0; // delete this pixel
				for (int i2=i-1; i2<=i+1; i2++) {
					for (int j2=j-1; j2<=j+1; j2++) {
						rm=( data[i2+j2*nx]==dt ? 1 : rm );
					}
				}
				if (rm>0){
					data2[i+jnx]=threshold;
					data[i+jnx]=dt+1;
				}
			}
		}
	}
	memcpy(data, data2, total_size * sizeof(float));
	if (verbose>0)
		printf("%d branches removed\n",nend);
	
	
	image->update();	
	if( data2 )
	{
		delete[]data2;
		data2 = 0;
	}
	
}

#ifdef SPARX_USING_CUDA

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
    std::cout<<"111111111  number of image==="<<N<<"number cluster=="<<K<<"img size"<<m<<std::endl;
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

    /*for (i = 0; i < K; i++) {
    	std::cout<<"average image"<<std::endl;
    	std::cout<<h_AVE[i*m+0]<<"  "<<h_AVE[i*m+1]<<"  "<<h_AVE[i*m+2]<<"  "<<h_AVE[i*m+3]<<"  "<<h_AVE[i*m+4]<<std::endl;
    }*/

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
// k-means SSE one iteration help function
/*int MPICUDA_kmeans::init_dist() {

    int status = cuda_mpi_kmeans_dist_SSE(h_AVE, d_AVE, h_dist, d_dist, d_im, h_im2, h_AVE2, h_asg, h_NC, params);

    return status;
}*/

/*int MPICUDA_kmeans::AVE_to_host() {

    int status = cuda_mpi_kmeans_copy_ave_from_device(h_AVE, d_AVE, params);

    return status;
}*/


int one_iter_SA();
// k-means SSE one iteration
/*int MPICUDA_kmeans::one_iter_SSE() {
    //if ( ite == 0)

    if( ite ==0) {
    	 int status_init=init_dist();
    	 ttt = compute_tt();//int status = 0;
    	 printf("intial energy ===%f\n",ttt);
    }
    //std::cout<<bb<<BB<<std::endl;
    int status = cuda_mpi_kmeans_SSE(h_AVE, d_AVE, h_dist, d_dist, d_im, h_im2, h_AVE2, h_asg, h_NC, params, ite, ttt);

    //float t = compute_tt();//int status = 0;
    //std::cout<<"engery at iteration"<<ite<<"==="<<t<<std::endl;

    ite++;
    return status;
}*/
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


float MPICUDA_kmeans::compute_tt() {
    /*int status = cuda_mpi_dist(h_AVE, d_AVE, h_dist, d_dist, d_im, n, K, m);
    vector <float> ji(K);
    int i;
    if (status != 0) {
        for (i = 0; i < K; ++i) ji[i] = -1.0;
        return -11111111;
    }
    for (i = 0; i < n; ++i) ji[h_asg[i]] += (h_im2[i] + h_AVE2[h_asg[i]] - 2 * h_dist[i * K + h_asg[i]]);
   float t =0.0;
   for (i =0; i<K;i++)
       t +=ji[i];
    return t;*/


    vector <float> ji(K);
    int i,j;
    float dist, temp;
    for (i = 0; i < n; ++i)
    {
    	dist =0;
	for ( j=0; j<m; j++)	 {
	    temp = (h_im[i*m+j] -h_AVE[ h_asg[i]*m+j]);
	    dist = dist + temp*temp;
	 }
    	ji[h_asg[i]] = ji[h_asg[i]]+ dist;
   }

   float t =0.0;
   for (i =0; i<K;i++)
       t +=ji[i];
    return t;
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

#endif //SPARX_USING_CUDA


