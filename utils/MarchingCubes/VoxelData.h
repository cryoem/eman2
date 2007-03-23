#ifndef _VOXEL_DATA_H_
#define _VOXEL_DATA_H_

#include "vecmath.h"
#include "reader.h"

class VoxelData {
public:
	VoxelData() : _init(false) {
		_volume = new Volume(2,2,2);
	}

	~VoxelData() {
		if(_init) {
			delete _volume;
			delete _reader;
		}
	}

	void loadMRC(const char* filename) {
		try {
			if(_init) {
				delete _volume;
				delete _reader;
			}
			_reader = new MRCReader(filename);
			_volume = _reader->getVolume();
			_volume->normalize(0,1000);
			_init = true;
		} catch(...) { _init = false; }
	}

	int getResolution() const {
		int resolution = 0;
		int num = 1;
		while(num < getSize()) {
			resolution++;
			num = 1 << resolution;
		}

		return resolution;
	}
	
	int getSize() const {
		return _volume->getSizeX();
	}

	float getPointValue(float x, float y, float z) {
		double fResult = 0.0;
        double fDx, fDy, fDz;
		Point3 source[] = {*(new Point3(0.5,0.5,0.5)), *(new Point3(0.8,0.1,0.1)), *(new Point3(0,0,0)), *(new Point3(0.2,0.8,0.8)), *(new Point3(0.5,0.1,0.1))};
		float weight[] = {0.5, 1.0, 1.5, 1.0, 0.5};
		for(int i=0; i<3; i++) {
			fDx = x - source[i][0];
			fDy = y - source[i][1];
			fDz = z - source[i][2];
			fResult += weight[i]/(fDx*fDx + fDy*fDy + fDz*fDz);
		}

        return fResult;
	}

	float getValue(int x, int y, int z) {
		//float val = _volume->getInterpDataAt(x*(_volume->getSizeX()-1), y*(_volume->getSizeY()-1), z*(_volume->getSizeZ()-1));
		//return _volume->getDataAt(x*(_volume->getSizeX()-1), y*(_volume->getSizeY()-1), z*(_volume->getSizeZ()-1));
		return _volume->getDataAt(x, y, z);
	}

	void setValue(float x, float y, float z, float value);

private:
	bool _init;
	Volume* _volume;
	VolumeReader* _reader;
};
#endif
