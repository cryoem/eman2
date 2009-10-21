// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Sasakthi S. Abeysinghe (sasakthi@gmail.com)
// Description:   Stores information of a density volume

#include "volume_data.h"

using namespace wustl_mm::SkeletonMaker;

		VolumeData::VolumeData(EMData* em) //eman2
		{
			this->emdata = em;
			SetSize( emdata->get_attr("nx"), emdata->get_attr("ny"), emdata->get_attr("nz") );
			SetSpacing( emdata->get_attr("apix_x"), emdata->get_attr("apix_y"), emdata->get_attr("apix_z") );
			SetOrigin( emdata->get_attr("origin_x"), emdata->get_attr("origin_y"), emdata->get_attr("origin_z") );
			owns_emdata = false;
		}

		VolumeData::VolumeData(int sizeX, int sizeY, int sizeZ) {
			InitializeVolumeData(sizeX, sizeY, sizeZ, 1, 1, 1, 0, 0, 0, true, 0);
		}

		VolumeData::VolumeData(int sizeX, int sizeY, int sizeZ, float val) {
			InitializeVolumeData(sizeX, sizeY, sizeZ, 1, 1, 1, 0, 0, 0, true, val);
		}

		VolumeData::VolumeData(int sizeX, int sizeY, int sizeZ, int offsetX, int offsetY, int offsetZ, VolumeData * data) {
			InitializeVolumeData(sizeX, sizeY, sizeZ, data->GetSpacingX(), data->GetSpacingY(), data->GetSpacingZ(), data->GetOriginX(), data->GetOriginY(), data->GetOriginZ(), false, 0);
			//int ct = 0 ;
			float value;//eman2
			for ( int k = offsetZ; k < sizeZ + offsetZ; k++) { //z in outer loop in eman2
				for (int j = offsetY; j < sizeY + offsetY; j++ ) {
					for (int i = offsetX; i < sizeX + offsetX; i++) {
						//this->data[ct] = data->GetDataAt(i, j, k);
						value = data->GetDataAt(i,j,k);//eman2
						SetDataAt(i-offsetX,j-offsetY,k-offsetZ, value);
						//ct++;
					}
				}
			}
		}

		VolumeData::VolumeData(VolumeData& obj)
		{
			size[0] = obj.GetSizeX();
			size[1] = obj.GetSizeY();
			size[2] = obj.GetSizeZ();
			spacing[0] = obj.GetSpacingX();
			spacing[1] = obj.GetSpacingY();
			spacing[2] = obj.GetSpacingZ();
			origin[0] = obj.GetOriginX();
			origin[1] = obj.GetOriginY();
			origin[2] = obj.GetOriginZ();
			emdata = new EMData( *obj.get_emdata() );
			owns_emdata = true;
		}
		VolumeData::~VolumeData() {
			cout << "~VolumeData()" << endl;
			//delete [] data;
			if (owns_emdata) //eman2
				delete emdata;
		}

		void VolumeData::InitializeVolumeData(int sizeX, int sizeY, int sizeZ, float spacingX, float spacingY, float spacingZ, float originX, float originY, float originZ, bool initializeData, float val) {
			emdata = new EMData(sizeX, sizeY, sizeZ); //eman2
			owns_emdata = true;
			SetSize(sizeX, sizeY, sizeZ);
			SetSpacing(spacingX, spacingY, spacingZ);
			SetOrigin(originX, originY, originZ);
			//int maxIndex = GetMaxIndex();
			//data = new float [maxIndex];
			if(initializeData) {
				emdata->to_value(val);
				//for(unsigned int i=0; i < maxIndex; i++) {
					//data[i] = val;
				//}
			}
		}

		int VolumeData::GetSize(int dimension) {
			return size[dimension];
		}

		int VolumeData::GetSizeX() {
			return GetSize(0);
		}

		int VolumeData::GetSizeY() {
			return GetSize(1);
		}

		int VolumeData::GetSizeZ() {
			return GetSize(2);
		}

		float VolumeData::GetSpacing(int dimension) {
			return spacing[dimension];
		}

		float VolumeData::GetSpacingX() {
			return GetSpacing(0);
		}

		float VolumeData::GetSpacingY() {
			return GetSpacing(1);
		}

		float VolumeData::GetSpacingZ() {
			return GetSpacing(2);
		}

		float VolumeData::GetOrigin(int dimension) {
			return origin[dimension];
		}

		float VolumeData::GetOriginX() {
			return GetOrigin(0);
		}

		float VolumeData::GetOriginY() {
			return GetOrigin(1);
		}

		float VolumeData::GetOriginZ() {
			return GetOrigin(2);
		}


		float VolumeData::GetDataAt(int x, int y, int z) {
			//return GetDataAt(GetIndex(x, y, z));
			return this->emdata->get_value_at(x, y, z);
		}

		float VolumeData::GetDataAt(int index) {
			//return data[index];
			return emdata->get_value_at(index);
		}

		int VolumeData::GetIndex(int x, int y, int z) {
			//return (x * GetSizeY() * GetSizeZ() + y * GetSizeZ() + z);
			return (x + y * GetSizeX() + z * GetSizeX() * GetSizeY());
		}

		int VolumeData::GetMaxIndex() {
			return size[0] * size[1] * size[2];
		}

		EMData * VolumeData::get_emdata() //eman2
		{
			return emdata;
		}

		void VolumeData::SetSpacing(float spacingX, float spacingY, float spacingZ) {
			// eman2 -->
			emdata->set_attr("apix_x", spacingX);
			emdata->set_attr("apix_y", spacingY);
			emdata->set_attr("apix_z", spacingZ);
			// <-- eman2
			spacing[0] = spacingX;
			spacing[1] = spacingY;
			spacing[2] = spacingZ;
		}

		void VolumeData::SetOrigin(float originX, float originY, float originZ) {
			// eman2 -->
			emdata->set_attr("origin_x", originX);
			emdata->set_attr("origin_y", originY);
			emdata->set_attr("origin_z", originZ);
			// <-- eman2
			origin[0] = originX;
			origin[1] = originY;
			origin[2] = originZ;
		}


		void VolumeData::SetSize(int sizeX, int sizeY, int sizeZ) {
			emdata->set_size(sizeX, sizeY, sizeZ); //eman2

			size[0] = sizeX;
			size[1] = sizeY;
			size[2] = sizeZ;
		}

		void VolumeData::SetDataAt(int x, int y, int z, float value) {
			//SetDataAt(GetIndex(x, y, z), value);
			emdata->set_value_at(x, y, z, value); //eman2
		}

		void VolumeData::SetDataAt(int index, float value) {
			//data[index] = value;
			//TODO: This may be a problem because EMAN and Gorgon do indexing differently --> see if this causes problems
			emdata->get_data()[index] = value; //eman2
		}
		void VolumeData::Pad(int padBy, double padValue) {

			int sizex = GetSizeX();
			int sizey = GetSizeY();
			int sizez = GetSizeZ();
			int newSizeX = sizex + 2*padBy;
			int newSizeY = sizey + 2*padBy;
			int newSizeZ = sizez + 2*padBy;

			// Method 1: Using get_clip ===========================================================
//			Region reg(0, 0, 0, newSizeX, newSizeY, newSizeZ);
//			EMData * new_emdata = emdata->get_clip(reg); //Extends the area x > sizex, y > sizey, z > sizez but does not translate (0,0,0)
//			new_emdata->translate(padBy,padBy,padBy); //We want equal padding on each side
//
//			for(int z = 0; z < newSizeZ; z++)
//				for(int y = 0; y < newSizeY; y++)
//					for(int x = 0; x < newSizeX; x++)
//						if ((x < padBy) || (y < padBy) || (z < padBy) || (x >= padBy + sizex) || (y >= padBy + sizey) || (z >= padBy + sizez))
//							new_emdata->set_value_at(x, y, z, padValue);
//
//			if (this->owns_emdata)
//				this->emdata->free_memory();
//			this->emdata = new_emdata;

			// Method 2: making a copy =============================================================
			float * newData = (float*) malloc ( newSizeX * newSizeY * newSizeZ * sizeof(float) );
			double value;

			for(int z = 0; z < newSizeZ; z++) {
				for(int y = 0; y < newSizeY; y++) {
					for(int x = 0; x < newSizeX; x++) {
						if ((x < padBy) || (y < padBy) || (z < padBy) || (x >= padBy + sizex) || (y >= padBy + sizey) || (z >= padBy + sizez)) {
							value = padValue;
						} else {
							value = GetDataAt(x-padBy, y-padBy, z-padBy);
						}

						//newData[x * newSizeY * newSizeZ + y * newSizeZ + z] = (float)value;
						newData[x + y * newSizeX + z * newSizeX * newSizeY] = static_cast<float>(value); //eman2
					}
				}
			}

			//code below changed for eman2
			float * data = emdata->get_data();
			SetSize(newSizeX, newSizeY, newSizeZ);
			emdata->set_data(newData, newSizeX, newSizeY, newSizeZ);
			//free(data); //I think this is handled by set_data

			// Method 3: doing in-place resize =====================================================
//			emdata->set_size(newSizeX, newSizeY, newSizeZ);
//			double val;
//			for (int k = newSizeZ - 1; k >= 0; k--)
//				for (int j = newSizeY - 1; k >= 0; j--)
//					for (int i = newSizeX - 1; k >=0; k--)
//					{
//						if ( i < padBy || i >= sizex+padBy || j < padBy || j >= sizey+padBy || k < padBy || k >= sizez+padBy)
//							emdata->set_value_at(i,j,k, padValue);
//						else
//						{
//							val = emdata->get_value_at(i-padBy, j-padBy, k-padBy);
//							emdata->set_value_at(i,j,k, float(val));
//						}
//					}
		}
