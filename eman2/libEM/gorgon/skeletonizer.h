// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Sasakthi S. Abeysinghe (sasakthi@gmail.com)
// Description:   Performs skeletonization on a grayscale volume

#include "volume.h"
using namespace wustl_mm::SkeletonMaker;

#ifndef GRAYSKELETONCPP_VOLUME_SKELETONIZER_H
#define GRAYSKELETONCPP_VOLUME_SKELETONIZER_H

const int MAX_GAUSSIAN_FILTER_RADIUS = 10;
const int DEFAULT_SKELETON_DIRECTION_RADIUS = 3;

namespace wustl_mm {
	namespace GraySkeletonCPP {

		class VolumeSkeletonizer
		{
		public:
			VolumeSkeletonizer(int pointRadius, int curveRadius, int surfaceRadius, int skeletonDirectionRadius=DEFAULT_SKELETON_DIRECTION_RADIUS);
			~VolumeSkeletonizer();
			static Volume * PerformPureJuSkeletonization(Volume * imageVol, string outputPath, double threshold, int minCurveWidth, int minSurfaceWidth);
			//Volume * PerformImmersionSkeletonizationAndPruning(Volume * sourceVol, Volume * preserveVol, double startGray, double endGray, double stepSize, int smoothingIterations, int smoothingRadius, int minCurveSize, int minSurfaceSize, int maxCurveHole, int maxSurfaceHole, string outputPath, bool doPruning, double pointThreshold, double curveThreshold, double surfaceThreshold);
			static void CleanUpSkeleton(Volume * skeleton, int minNumVoxels = 4, float valueThreshold = 0.5); //Added for EMAN2
			static void MarkSurfaces(Volume* skeleton); //Added for EMAN2
			
		private:
			static bool Are26Neighbors(Vec3<int> u, Vec3<int> v); //Added for EMAN2
			static Volume * GetJuSurfaceSkeleton(Volume * sourceVolume, Volume * preserve, double threshold);
			static Volume * GetJuCurveSkeleton(Volume * sourceVolume, Volume * preserve, double threshold, bool is3D);
			static Volume * GetJuTopologySkeleton(Volume * sourceVolume, Volume * preserve, double threshold);
			static void PruneCurves(Volume * sourceVolume, int pruneLength);
			static void PruneSurfaces(Volume * sourceVolume, int pruneLength);
			static void VoxelOr(Volume * sourceAndDestVolume1, Volume * sourceVolume2);
			static Volume * GetJuThinning(Volume * sourceVolume, Volume * preserve, double threshold, char thinningClass);

			static const char THINNING_CLASS_SURFACE_PRESERVATION;
			static const char THINNING_CLASS_CURVE_PRESERVATION_2D;
			static const char THINNING_CLASS_CURVE_PRESERVATION;
			static const char THINNING_CLASS_POINT_PRESERVATION;
			static const char THINNING_CLASS_TOPOLOGY_PRESERVATION;
			static const char PRUNING_CLASS_PRUNE_SURFACES;
			static const char PRUNING_CLASS_PRUNE_CURVES;
			static const char PRUNING_CLASS_PRUNE_POINTS;

			//MathLib * math;
			//NormalFinder * surfaceNormalFinder;
			//ProbabilityDistribution3D gaussianFilterPointRadius;
			//ProbabilityDistribution3D gaussianFilterCurveRadius;
			//ProbabilityDistribution3D gaussianFilterSurfaceRadius;
			//ProbabilityDistribution3D gaussianFilterMaxRadius;
			//ProbabilityDistribution3D uniformFilterSkeletonDirectionRadius;
			int pointRadius;
			int curveRadius;
			int surfaceRadius;
			int skeletonDirectionRadius;
		};
	}
}

#endif
