#ifndef _MESH_H_
#define _MESH_H_

#include "vecmath.h"

#include <vector>
//#include <hash_map>
#include <map>

#include <GL/gl.h>
#include <GL/glu.h>
//#include <GL/glut.h>

using std::map;
//using std::hash_map;
using std::vector;

// /SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup
namespace EMAN
{

	class Mesh {
	public:
	
		Mesh(vector<float>* points, vector<int>* faces, vector<float>* normals) : _points(points), _faces(faces), _normals(normals) {
			_normalsSm = new vector<float>();
			// Compute smooth normals.
			map<int, vector<int> > point_map;
	
			for(unsigned int i=0; i<_faces->size(); i+=3) {
				for(int j=0; j<3; j++) {
					int pt = (*_faces)[i+j];
					if(!point_map.count(pt)) 
						point_map[pt] = vector<int>();
					point_map[pt].push_back(i);
				}
			}
	
			for(unsigned int i=0; i<_points->size(); i+=3) {
				int size = point_map[i].size();
				float x,y,z;
				x = y = z = 0;
				Vector3 n(0,0,0);
				//std::cout << size << std::endl;
				for(int j=0; j<size; j++) {
					int index = point_map[i][j];
					Vector3 comp((*_normals)[index], (*_normals)[index+1], (*_normals)[index+2]);
					n += comp;
					//x += (*_normals)[index];
					//y += (*_normals)[index+1];
					//z += (*_normals)[index+2];
				}	
				n.normalize();
				//_normalsSm->push_back(x/((float)size));
				//_normalsSm->push_back(y/((float)size));
				//_normalsSm->push_back(z/((float)size));
				_normalsSm->push_back(n[0]);
				_normalsSm->push_back(n[1]);
				_normalsSm->push_back(n[2]);
			}
		}
	
		~Mesh() {
			delete _points;
			delete _faces;
			delete _normals;
			delete _normalsSm;
		}
	
	
		void draw(bool smooth = false) const {
			glBegin(GL_TRIANGLES);
			int size = _faces->size();
			for(int i=0; i<size; i+=3) {
				if(!smooth)
					glNormal3f((*_normals)[i], (*_normals)[i+1], (*_normals)[i+2]);
				for(int j=0; j<3; j++) {
					int pt = (*_faces)[i+j];
					if(smooth)
						glNormal3f((*_normalsSm)[pt], (*_normalsSm)[pt+1], (*_normalsSm)[pt+2]);
	
					glVertex3f((*_points)[pt], (*_points)[pt+1], (*_points)[pt+2]);
				}
			}
			glEnd();
			
		}
	
	private:
		vector<float>* _points;
		vector<int>* _faces;
		vector<float>* _normals;
		vector<float>* _normalsSm;
	};

}

#endif
