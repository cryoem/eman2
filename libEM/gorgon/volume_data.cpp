// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Sasakthi S. Abeysinghe (sasakthi@gmail.com)
// Description:   Stores information of a density volume

#include "volume_data.h"

using namespace wustl_mm::SkeletonMaker;

		VolumeData::VolumeData(EMData* em)
		{
			this->emdata = em;
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
			float value;
			for ( int k = offsetZ; k < sizeZ + offsetZ; k++) {
				for (int j = offsetY; j < sizeY + offsetY; j++ ) {
					for (int i = offsetX; i < sizeX + offsetX; i++) {
						value = data->GetDataAt(i,j,k);
						SetDataAt(i-offsetX,j-offsetY,k-offsetZ, value);
					}
				}
			}
		}

		VolumeData::VolumeData(VolumeData& obj)
		{
			emdata = new EMData( *obj.get_emdata() );
			owns_emdata = true;
			//SetSize( obj.GetSizeX(), obj.GetSizeY(), obj.GetSizeZ() );
			//SetSpacing( obj.GetSpacingX(), obj.GetSpacingY(), obj.GetSpacingZ() );
			//SetOrigin( obj.GetOriginX(), GetOriginY(), obj.GetOriginZ() );
		}
		VolumeData::~VolumeData() {
			if (owns_emdata)
				delete emdata;
		}

		void VolumeData::InitializeVolumeData(int sizeX, int sizeY, int sizeZ, float spacingX, float spacingY, float spacingZ, float originX, float originY, float originZ, bool initializeData, float val) {
			emdata = new EMData(sizeX, sizeY, sizeZ);
			owns_emdata = true;
			//SetSize(sizeX, sizeY, sizeZ);
			SetSpacing(spacingX, spacingY, spacingZ);
			SetOrigin(originX, originY, originZ);
			if(initializeData) {
				emdata->to_value(val);
			}
		}

		int VolumeData::GetSize(int dimension) {
			int ret = 0;
			switch (dimension)
			{
				case 0: ret = GetSizeX();
					break;
				case 1: ret = GetSizeY();
					break;
				case 2: ret = GetSizeZ();
					break;
				default:
					throw InvalidParameterException("VolumeData::GetSize requires an argument of 0, 1, or 2");
			}

			return ret;

		}

		int VolumeData::GetSizeX() {
			return emdata->get_xsize();
		}

		int VolumeData::GetSizeY() {
			return emdata->get_ysize();
		}

		int VolumeData::GetSizeZ() {
			return emdata->get_zsize();
		}

		float VolumeData::GetSpacing(int dimension) {
			float ret = 0;
			switch (dimension)
			{
				case 0: ret = GetSpacingX();
					break;
				case 1: ret = GetSpacingY();
					break;
				case 2: ret = GetSpacingZ();
					break;
				default:
					throw InvalidParameterException("VolumeData::GetSpacing requires an argument of 0, 1, or 2");
			}

			return ret;
		}

		float VolumeData::GetSpacingX() {
			return emdata->get_attr("apix_x");
		}

		float VolumeData::GetSpacingY() {
			return emdata->get_attr("apix_y");
		}

		float VolumeData::GetSpacingZ() {
			return emdata->get_attr("apix_y");
		}

		float VolumeData::GetOrigin(int dimension) {
			float ret = 0;
			switch (dimension)
			{
				case 0: ret = GetOriginX();
					break;
				case 1: ret = GetOriginY();
					break;
				case 2: ret = GetOriginZ();
					break;
				default:
					throw InvalidParameterException("VolumeData::GetOrigin requires an argument of 0, 1, or 2");
			}

			return ret;
		}

		float VolumeData::GetOriginX() {
			return emdata->get_attr("origin_x");
		}

		float VolumeData::GetOriginY() {
			return emdata->get_attr("origin_y");
		}

		float VolumeData::GetOriginZ() {
			return emdata->get_attr("origin_z");
		}


		float VolumeData::GetDataAt(int x, int y, int z) {
			return this->emdata->get_value_at(x, y, z);
		}

		float VolumeData::GetDataAt(int index) {
			//TODO: This may be a problem because EMAN2 and Gorgon do indexing differently --> see if this causes problems
			return emdata->get_value_at(index);
		}

		int VolumeData::GetIndex(int x, int y, int z) {
			return (x + y * GetSizeX() + z * GetSizeX() * GetSizeY());
		}

		int VolumeData::GetMaxIndex() {
			return GetSizeX() * GetSizeY() * GetSizeZ();
		}

		EMData * VolumeData::get_emdata() //eman2
		{
			return emdata;
		}

		void VolumeData::SetSpacing(float spacingX, float spacingY, float spacingZ) {
			emdata->set_attr("apix_x", spacingX);
			emdata->set_attr("apix_y", spacingY);
			emdata->set_attr("apix_z", spacingZ);
		}

		void VolumeData::SetOrigin(float originX, float originY, float originZ) {
			emdata->set_attr("origin_x", originX);
			emdata->set_attr("origin_y", originY);
			emdata->set_attr("origin_z", originZ);
		}


		void VolumeData::SetSize(int sizeX, int sizeY, int sizeZ) {
			emdata->set_size(sizeX, sizeY, sizeZ);
		}

		void VolumeData::SetDataAt(int x, int y, int z, float value) {
			emdata->set_value_at(x, y, z, value);
		}

		void VolumeData::SetDataAt(int index, float value) {
			//TODO: This may be a problem because EMAN2 and Gorgon do indexing differently --> see if this causes problems
			emdata->get_data()[index] = value;
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

						newData[x + y * newSizeX + z * newSizeX * newSizeY] = static_cast<float>(value); //eman2
					}
				}
			}

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
