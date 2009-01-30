#ifndef cuda_project_h__
#define cuda_project_h__ 1


#include "transform.h"
using EMAN::Transform;

float** main_t(const float* ,const int,const int,const int);

float* standard_project(const Transform* const t, const int nx, const int ny);

#endif //  cuda_project_h__ 1