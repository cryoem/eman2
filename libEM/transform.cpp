/**
 * $Id$
 */
#include "transform.h"
#include "exception.h"

using namespace EMAN;

const float Transform::ERR_LIMIT = 0.000001f;


Transform::Transform()
{
	init();
	to_identity();
}



Transform::Transform(const vector<float>& rotation, EulerType euler_type)
{
	init();
	set_rotation(rotation, euler_type);
}


Transform::Transform(const vector<float>& rotation, EulerType euler_type, 
					 const Vec3f& posttrans)
{
	init();
	set_rotation(rotation, euler_type);
	set_posttrans(posttrans);
}

		
Transform::Transform(const vector<float>& rotation, EulerType euler_type, 
					 const Vec3f & pretrans, const Vec3f& posttrans)
{
	init();
	set_rotation(rotation, euler_type);
	set_pretrans(pretrans);
	set_posttrans(posttrans);
}

Transform::Transform(map<string, float>& rotation, EulerType euler_type)
{
	init();
	set_rotation(rotation, euler_type);
}


Transform::Transform(map<string, float>& rotation, EulerType euler_type, 
					 const Vec3f& posttrans)
{
	init();
	set_rotation(rotation, euler_type);
	set_posttrans(posttrans);
}

		
Transform::Transform(map<string, float>& rotation, EulerType euler_type, 
					 const Vec3f & pretrans, const Vec3f& posttrans)
{
	init();
	set_rotation(rotation, euler_type);
	set_pretrans(pretrans);
	set_posttrans(posttrans);
}

		
Transform::~Transform()
{
	const int n = 4;
	for (int i = 0; i < n; i++) {
		delete [] matrix[i];
		matrix[i] = 0;
	}
	delete [] matrix;
}


void Transform::init()
{
	const int n = 4;
	const int m = 3;
	
	matrix = new float*[4];
	for (int i = 0; i < n; i++) {
		matrix[i] = new float[m];
		memset(matrix[i], 0, sizeof(float) * m);
	}
}


void Transform::to_identity()
{
	for (int i = 0; i < 3; i++) {
		matrix[i][i] = 1;
	}	
}


void Transform::orthogonalize()	
{

}

Transform Transform::inverse()
{
	return Transform();
}

	
void Transform::set_pretrans(const Vec3f & pretrans)
{
	pre_trans = pretrans;
}

void Transform::set_posttrans(const Vec3f & posttrans)
{
	for (int i = 0; i < 3; i++) {
		matrix[3][i] = posttrans[i];
	}
}

void Transform::set_center(const Vec3f & center)
{
	pre_trans= -center;
	for (int i = 0; i < 3; i++) {
		matrix[3][i]=center[i];
	}
}

void Transform::set_rotation(const vector<float> & rotation, EulerType euler_type)
{
	float a0 = rotation[0];
	float a1 = rotation[1];
	float a2 = rotation[2];
	float a3 = 0;
	if (rotation.size() > 3) {
		a3 = rotation[3];
	}
	
	bool is_quaternion = 0;

	switch(euler_type) {
	case EMAN:
	case IMAGIC:
		break;
		
	case SPIDER:
		a1 = a1 - 1.5f * M_PI;
		a2 = a2 - 0.5f * M_PI;
		break;
		
	case MRC:
		a1 = fmod(-a1 + 2.5f * M_PI, 2 * M_PI);
		a2 = fmod(a2 + 0.5f * M_PI, 2.0f * M_PI);
		break;
		
	case QUATERNION:
		is_quaternion = 1;
		break;
		
	case SPIN:
	case SGIROT:
		{
			is_quaternion = 1;
			float f = sin(a0 / 2.0f);
			a1 *= f;
			a2 *= f;
			a3 *= f;
			a0 = cos(a0 / 2.0f);
		}
		break;
		
	case MATRIX:
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				matrix[i][j] = rotation[i*3+j];
			}
		}
		break;
		
	default:
		throw InvalidValueException(euler_type, "unknown Euler Type");
	}
	
	if (is_quaternion) {
		quaternion2matrix(a0, a1, a2, a3);
	}
	else if (euler_type != MATRIX) {
		matrix[0][0] = cos(a3)*cos(a2) - cos(a1)*sin(a2)*sin(a3);
		matrix[0][1] = -(sin(a3)*cos(a2) + cos(a1)*sin(a2)*cos(a3));
		matrix[0][2] = sin(a1)*sin(a2);
		matrix[1][0] = cos(a3)*sin(a2) + cos(a1)*cos(a2)*sin(a3);
		matrix[1][1] = -sin(a3)*sin(a2) + cos(a1)*cos(a2)*cos(a3);
		matrix[1][2] = -sin(a1)*cos(a2);
		matrix[2][0] = sin(a1)*sin(a3);
		matrix[2][1] = sin(a1)*cos(a3);
		matrix[2][2] = cos(a1);
	}
	
}

void Transform::set_rotation(map<string, float> & rotation, EulerType euler_type)
{
	float a0 = 0;
	float a1 = 0;
	float a2 = 0;
	float a3 = 0;
	bool is_quaternion = 0;

	switch(euler_type) {
	case EMAN:
		a1 = rotation["alt"];
		a2 = rotation["az"];
		a3 = rotation["phi"];
		break;
	case IMAGIC:
		a1 = rotation["alpha"];
		a2 = rotation["beta"];
		a3 = rotation["gamma"];
		break;
		
	case SPIDER:
		a1 = rotation["phi"];
		a2 = rotation["theta"] - 1.5f * M_PI;
		a3 = rotation["gamma"] - 0.5f * M_PI;
		break;
		
	case MRC:
		a1 = rotation["theta"];
		a2 = fmod(-rotation["phi"] + 2.5f * M_PI, 2.0 * M_PI);
		a3 = fmod(rotation["omega"] + 0.5f * M_PI, 2.0f * M_PI);
		break;
		
	case QUATERNION:
		a0 = rotation["e0"];
		a1 = rotation["e1"];
		a2 = rotation["e2"];
		a3 = rotation["e3"];
		
		is_quaternion = 1;
		break;
		
	case SPIN:
	case SGIROT:
		{
			a0 = rotation["q"];
			a1 = rotation["n1"];
			a2 = rotation["n2"];
			a3 = rotation["n3"];
		
			is_quaternion = 1;
			float f = sin(a0 / 2.0f);
			a1 *= f;
			a2 *= f;
			a3 *= f;
			a0 = cos(a0 / 2.0f);
		}
		break;
		
	case MATRIX:
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				char str[32];
				sprintf(str, "m%0d%0d", i, j);
				matrix[i][j] = rotation[str];
			}
		}
		break;
		
	default:
		throw InvalidValueException(euler_type, "unknown Euler Type");
	}
	
	if (is_quaternion) {
		quaternion2matrix(a0, a1, a2, a3);
	}
	else if (euler_type != MATRIX) {
		matrix[0][0] =  cos(a3)*cos(a2) - cos(a1)*sin(a2)*sin(a3);
		matrix[0][1] = -sin(a3)*cos(a2) - cos(a1)*sin(a2)*cos(a3);
		matrix[0][2] =  sin(a1)*sin(a2);
		matrix[1][0] =  cos(a3)*sin(a2) + cos(a1)*cos(a2)*sin(a3);
		matrix[1][1] = -sin(a3)*sin(a2) + cos(a1)*cos(a2)*cos(a3);
		matrix[1][2] = -sin(a1)*cos(a2);
		matrix[2][0] =  sin(a1)*sin(a3);
		matrix[2][1] =  sin(a1)*cos(a3);
		matrix[2][2] =  cos(a1);
	}
	
}

void Transform::set_scale(float scale)
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 3; j++) {
			matrix[i][j] *= scale;
		}
	}
}

Vec3f Transform::get_center() const
{
	return Vec3f();
}

Vec3f Transform::get_pretrans() const
{
	return pre_trans;
}

Vec3f Transform::get_posttrans() const
{
	return Vec3f(matrix[3][0], matrix[3][1], matrix[3][2]);
}


float Transform::eman_alt() const
{
	float alt = 0;
	float max = 1 - ERR_LIMIT;
	float sca=get_scale();
	float mx=matrix[2][2]/sca;
	
	if (mx > max) {
		alt = 0;
	}
	else if (mx < -max) {
		alt = M_PI;
	}
	else {
		alt = (float) acos(mx);
	}
	return alt;
}

float Transform::eman_az() const
{
	float az = 0;
	float max = 1 - ERR_LIMIT;
	float sca=get_scale();
	float mx=matrix[2][2]/sca;
	if (fabs(mx) > max) {
		az = (float)(atan2(matrix[1][0], matrix[0][0]);
	}
	else {
		az = (float)atan2(matrix[0][2], -matrix[1][2]);
	}
	return az;
}

float Transform::eman_phi() const
{
	float phi = 0;
	float sca=get_scale();

	if (fabs(matrix[2][2]/sca) > (1 - ERR_LIMIT)) {
		phi = 0;
	}
	else {
		phi = (float)atan2(matrix[2][0], matrix[2][1]);
	}

	return phi;
}

	

map<string,float> Transform::get_rotation(EulerType euler_type) const
{
	map<string, float> result;
	
	float alt = eman_alt();
	float az = eman_az();
	float phi = eman_phi();
	
	switch (euler_type) {
	case EMAN:
		result["alt"] = alt;
		result["az"] = az;
		result["phi"] = phi;
		break;

	case MRC:
		result["theta"] = alt;
		result["phi"] = fmod(-az + 2.5f * M_PI, 2 * M_PI);
		result["omega"] = fmod(phi + 1.5f * M_PI, 2 * M_PI);
		break;

	case IMAGIC:
		result["alpha"] = phi;
		result["beta"] = alt;
		result["gamma"] = az;
		break;

	case SPIDER:
		result["phi"] = fmod((az + 1.5f * M_PI), 2 * M_PI);
		result["theta"] = alt;
		result["gamma"] = fmod((phi + M_PI/2), 2 * M_PI);	
		break;

	case SPIN:
	case SGIROT:
		{
			vector<float> sgivec = matrix2sgi();
			result["q"] = sgivec[0];
			result["n1"] = sgivec[1];
			result["n2"] = sgivec[2];
			result["n3"] = sgivec[3];
		}
		break;
		
	case QUATERNION:
		{
			vector<float> quatvec = matrix2quaternion();
			result["e0"] = quatvec[0];
			result["e1"] = quatvec[1];
			result["e2"] = quatvec[2];
			result["e3"] = quatvec[3];
		}
		break;
		
	case MATRIX:
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				char str[32];
				sprintf(str, "m%d%d", i, j);
				result[string(str)] = matrix[i][j];
			}
		}
		break;
	default:
		throw InvalidValueException(euler_type, "unknown Euler Type");
	}

	return result;
}

float Transform::get_scale() const
{
	// Assumes uniform scaling, calculation uses Z only
	float ksq=SQR(matrix[0][2])+SQR(matrix[1][2])+SQR(matrix[2][2]);	// should be 1.0 if no scaling

	return sqrt(ksq);
}

		
float Transform::orthogonality() const
{
	return 0;
}

		
float * Transform::operator[] (int i)
{
	return matrix[i];
}

const float * Transform::operator[] (int i) const
{
	return matrix[i];
}

vector<float> Transform::matrix2quaternion() const
{
	vector<float> result(4);
	int i = 0;

	if (matrix[0][0] > matrix[1][1]) {
		if (matrix[0][0] > matrix[2][2]) {
			i = 0;
		}
		else {
			i = 2;
		}
	}
	else {
		if (matrix[1][1] > matrix[2][2]) {
			i = 1;
		}
		else {
			i = 2;
		}
	}

	if (matrix[0][0] + matrix[1][1] + matrix[2][2] > matrix[i][i]) {
		result[0] = (float) (sqrt(matrix[0][0] + matrix[1][1] + matrix[2][2] + 1) / 2.0);
		result[1] = (float) ((matrix[1][2] - matrix[2][1]) / (4 * result[0]));
		result[2] = (float) ((matrix[2][0] - matrix[0][2]) / (4 * result[0]));
		result[3] = (float) ((matrix[0][1] - matrix[1][0]) / (4 * result[0]));
	}
	else {
		float quat[3];
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;

		quat[i] = (float) (sqrt(matrix[i][i] - matrix[j][j] - matrix[k][k] + 1) / 2.0);
		quat[j] = (float) ((matrix[i][j] + matrix[j][i]) / (4 * quat[i]));
		quat[k] = (float) ((matrix[i][k] + matrix[k][i]) / (4 * quat[i]));

		result[0] = (float) ((matrix[j][k] - matrix[k][j]) / (4 * quat[i]));
		result[1] = quat[0];
		result[2] = quat[1];
		result[3] = quat[2];
	}
							 
	return result;
}


vector<float> Transform::matrix2sgi() const
{
	vector<float> q = matrix2quaternion();
	Vec3f vec(q[1], q[2], q[3]);
	float len = vec.length();
	float radians = 0;
	Vec3f axis;
	
	if (len > 0.00001f) {
		radians = 2.0f * acos(q[0]);
		axis = vec * ((float) (1.0f / len));
	}
	else {
		radians = 0;
		axis.set_value(0.0f, 0.0f, 1.0f);
	}

	vector<float> result(4);
	result[0] = radians;
	result[1] = axis[0];
	result[2] = axis[1];
	result[3] = axis[2];
	
	return result;
}

void Transform::quaternion2matrix(float e0, float e1, float e2, float e3)
{
	matrix[0][0] = e0 * e0 + e1 * e1 - e2 * e2 - e3 * e3;
	matrix[0][1] = 2.0f * (e1 * e2 + e0 * e3);
	matrix[0][2] = 2.0f * (e1 * e3 - e0 * e2);

	matrix[1][0] = 2.0f * (e1 * e2 - e0 * e3);
	matrix[1][1] = e0 * e0 + e1 * e1 + e2 * e2 - e3 * e3;
	matrix[1][2] = 2.0f * (e2 * e3 + e0 * e1);

	matrix[2][0] = 2.0f * (e1 * e3 + e0 * e2);
	matrix[2][1] = 2.0f * (e2 * e3 - e0 * e1);
	matrix[2][2] = e0 * e0 - e1 * e1 - e2 * e2 + e3 * e3;
}

