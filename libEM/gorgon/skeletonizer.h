// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Sasakthi S. Abeysinghe (sasakthi@gmail.com)
// Description:   Performs skeletonization on a grayscale volume

#include "volume.h"
using namespace wustl_mm::SkeletonMaker;

#ifndef GRAYSKELETONCPP_VOLUME_SKELETONIZER_H
#define GRAYSKELETONCPP_VOLUME_SKELETONIZER_H

const int MAX_GAUSSIAN_FILTER_RADIUS = 10;

namespace wustl_mm {
	namespace GraySkeletonCPP {
		class VolumeSkeletonizer
		{
		public:
			Volume * PerformPureJuSkeletonization(Volume * imageVol, string outputPath, double threshold, int minCurveWidth, int minSurfaceWidth);
			//Volume * PerformImmersionSkeletonizationAndPruning(Volume * sourceVol, Volume * preserveVol, double startGray, double endGray, double stepSize, int smoothingIterations, int smoothingRadius, int minCurveSize, int minSurfaceSize, int maxCurveHole, int maxSurfaceHole, string outputPath, bool doPruning, double pointThreshold, double curveThreshold, double surfaceThreshold);
		private:
			Volume * GetJuSurfaceSkeleton(Volume * sourceVolume, Volume * preserve, double threshold);
			Volume * GetJuCurveSkeleton(Volume * sourceVolume, Volume * preserve, double threshold, bool is3D);
			Volume * GetJuTopologySkeleton(Volume * sourceVolume, Volume * preserve, double threshold);
			void PruneCurves(Volume * sourceVolume, int pruneLength);
			void PruneSurfaces(Volume * sourceVolume, int pruneLength);
			void VoxelOr(Volume * sourceAndDestVolume1, Volume * sourceVolume2);
			Volume * GetJuThinning(Volume * sourceVolume, Volume * preserve, double threshold, char thinningClass);

			static const char THINNING_CLASS_SURFACE_PRESERVATION;
			static const char THINNING_CLASS_CURVE_PRESERVATION_2D;
			static const char THINNING_CLASS_CURVE_PRESERVATION;
			static const char THINNING_CLASS_POINT_PRESERVATION;
			static const char THINNING_CLASS_TOPOLOGY_PRESERVATION;
			static const char PRUNING_CLASS_PRUNE_SURFACES;
			static const char PRUNING_CLASS_PRUNE_CURVES;
			static const char PRUNING_CLASS_PRUNE_POINTS;
		};
	}
}

#endif
