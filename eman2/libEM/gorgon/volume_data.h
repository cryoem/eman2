// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Sasakthi S. Abeysinghe (sasakthi@gmail.com)
// Description:   Stores information of a density volume
#include "emdata.h"
using namespace EMAN;

#ifndef SKELETON_MAKER_VOLUME_DATA_H
#define SKELETON_MAKER_VOLUME_DATA_H

namespace wustl_mm {
	namespace SkeletonMaker {
		class VolumeData {
		public:
		//Constructors & Destructor
			VolumeData(EMData* em); //eman2
			VolumeData(int sizeX, int sizeY, int sizeZ);
			VolumeData(int sizeX, int sizeY, int sizeZ, float val);
			VolumeData(int sizeX, int sizeY, int sizeZ, int offsetX, int offsetY, int offsetZ, VolumeData * data);
			VolumeData(VolumeData& obj); //eman2
			~VolumeData();

		//Member functions
			int GetSize(int dimension);
			int GetSizeX();
			int GetSizeY();
			int GetSizeZ();
			float GetSpacing(int dimension);
			float GetSpacingX();
			float GetSpacingY();
			float GetSpacingZ();
			float GetOrigin(int dimension);
			float GetOriginX();
			float GetOriginY();
			float GetOriginZ();
			float GetDataAt(int x, int y, int z);
			float GetDataAt(int index);
			int GetIndex(int x, int y, int z);
			int GetMaxIndex();

			void SetSpacing(float spacingX, float spacingY, float spacingZ);
			void SetOrigin(float originX, float originY, float originZ);
			void SetDataAt(int x, int y, int z, float value);
			void SetDataAt(int index, float value);
			void Pad(int padBy, double padValue);
			EMData * get_emdata(); //eman2
		private:
			void InitializeVolumeData(int sizeX, int sizeY, int sizeZ, float spacingX, float spacingY, float spacingZ, float originX, float originY, float originZ, bool initializeData, float val);
			void SetSize(int sizeX, int sizeY, int sizeZ);

		//Member variables
		public:
			bool owns_emdata; //eman2
		private:
			EMData * emdata; //eman2
			
		};


	}
}


#endif
