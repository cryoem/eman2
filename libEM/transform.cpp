/**
 * $Id$
 */
#include "transform.h"
#include "log.h"
#include <math.h>
#include <ctype.h>
#ifdef WIN32
#define M_PI 3.14159265358979323846
#endif
using namespace EMAN;

Matrix3f::Matrix3f(float m0, float m1, float m2, float m3, float m4, float m5, float m6, float m7,
				   float m8)
{
	matrix = gsl_matrix_alloc(3, 3);

	gsl_matrix_set(matrix, 0, 0, m0);
	gsl_matrix_set(matrix, 0, 1, m1);
	gsl_matrix_set(matrix, 0, 2, m2);

	gsl_matrix_set(matrix, 1, 0, m3);
	gsl_matrix_set(matrix, 1, 1, m4);
	gsl_matrix_set(matrix, 1, 2, m5);

	gsl_matrix_set(matrix, 2, 0, m6);
	gsl_matrix_set(matrix, 2, 1, m7);
	gsl_matrix_set(matrix, 2, 2, m8);

}

Matrix3f::Matrix3f(const vector < float >&m)
{
	const int n = 3;
	matrix = gsl_matrix_alloc(n, n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			gsl_matrix_set(matrix, i, j, m[i * n + j]);
		}
	}
}

Matrix3f::Matrix3f(gsl_matrix * m)
{
	const int n = 3;
	matrix = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(matrix, m);
}

Matrix3f::Matrix3f(const Matrix3f & m)
{
	const int n = 3;
	matrix = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(matrix, m.get_gsl_matrix());
}

Matrix3f::~Matrix3f()
{
	if (matrix) {
		gsl_matrix_free(matrix);
		matrix = 0;
	}
}


Matrix3f & Matrix3f::operator=(const Matrix3f & m)
{
	if (this != &m) {
		if (matrix) {
			gsl_matrix_free(matrix);
		}
		const int n = 3;
		matrix = gsl_matrix_alloc(n, n);
		gsl_matrix_memcpy(matrix, m.get_gsl_matrix());
	}
	return *this;
}

void Matrix3f::make_identity()
{
	gsl_matrix_set_identity(matrix);
}

Matrix3f & Matrix3f::mult_right(const Matrix3f & m)
{
	(*this) = (*this) * m;
	return *this;
}

Matrix3f & Matrix3f::mult_left(const Matrix3f & m)
{
	Matrix3f mcopy = m;
	mcopy = mcopy * (*this);
	(*this) = mcopy;
	return *this;
}

gsl_matrix *Matrix3f::get_gsl_matrix() const
{
	return matrix;
}

void Matrix3f::set_value(const vector < float >&m)
{
	const int n = 3;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			gsl_matrix_set(matrix, i, j, m[i * n + j]);
		}
	}
}

vector < float >Matrix3f::get_as_list() const
{
	const int n = 3;
	vector < float >m(n * n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			m[i * n + j] = gsl_matrix_get(matrix, i, j);
		}
	}
	return m;
}

// todo: check whether the matrix can be inversed or not.
Matrix3f & Matrix3f::inverse()
{
	const int n = 3;
	gsl_matrix *inverse = gsl_matrix_alloc(n, n);
	gsl_permutation *p = gsl_permutation_alloc(n);

	int s = 0;
	gsl_linalg_LU_decomp(matrix, p, &s);
	gsl_linalg_LU_invert(matrix, p, inverse);
	gsl_matrix_memcpy(matrix, inverse);
	gsl_matrix_free(inverse);
	inverse = 0;

	return (*this);
}

Matrix3f & Matrix3f::transpose()
{
	gsl_matrix_transpose(matrix);
	return (*this);
}

Matrix3f Matrix3f::create_inverse() const
{
	const int n = 3;
	gsl_matrix *matrix_copy = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(matrix_copy, matrix);

	gsl_matrix *inverse = gsl_matrix_alloc(n, n);
	gsl_permutation *p = gsl_permutation_alloc(n);

	int s = 0;
	gsl_linalg_LU_decomp(matrix_copy, p, &s);
	gsl_linalg_LU_invert(matrix_copy, p, inverse);

	gsl_matrix_free(matrix_copy);
	matrix_copy = 0;

	return Matrix3f(inverse);
}

Vec3 < float >Matrix3f::get_row(int i) const
{
	return Vec3 < float >(gsl_matrix_get(matrix, i, 0),
						  gsl_matrix_get(matrix, i, 1),
						  gsl_matrix_get(matrix, i, 2));

}


Vec3 < float >Matrix3f::get_col(int i) const
{
	return Vec3 < float >(gsl_matrix_get(matrix, 0, i),
						  gsl_matrix_get(matrix, 1, i),
						  gsl_matrix_get(matrix, 2, i));

}


Matrix3f EMAN::operator+(float f, const Matrix3f & m2)
{
	Matrix3f m = m2;
	m += f;
	return m;
}

Matrix3f EMAN::operator-(float f, const Matrix3f & m2)
{
	Matrix3f m = m2;
	m -= f;
	return m;
}

Matrix3f EMAN::operator*(float scale, const Matrix3f & m2)
{
	Matrix3f m = m2;
	m *= scale;
	return m;
}

Matrix3f EMAN::operator/(float scale, const Matrix3f & m2)
{
	Matrix3f m = m2;
	m /= scale;
	return m;
}


Matrix3f EMAN::operator+(const Matrix3f & m1, float f)
{
	Matrix3f m = m1;
	m += f;
	return m;
}

Matrix3f EMAN::operator-(const Matrix3f & m1, float f)
{
	Matrix3f m = m1;
	m -= f;
	return m;
}

Matrix3f EMAN::operator*(const Matrix3f & m1, float scale)
{
	Matrix3f m = m1;
	m *= scale;
	return m;
}

Matrix3f EMAN::operator/(const Matrix3f & m1, float scale)
{
	Matrix3f m = m1;
	m /= scale;
	return m;
}


Matrix3f EMAN::operator+(const Matrix3f & m1, const Matrix3f & m2)
{
	Matrix3f m = m1;
	m += m2;
	return m;
}

Matrix3f EMAN::operator-(const Matrix3f & m1, const Matrix3f & m2)
{
	Matrix3f m = m1;
	m -= m2;
	return m;
}

Matrix3f EMAN::operator*(const Matrix3f & m1, const Matrix3f & m2)
{
	Matrix3f m = m1;
	m *= m2;
	return m;
}

Matrix3f EMAN::operator/(const Matrix3f & m1, const Matrix3f & m2)
{
	Matrix3f m = m1;
	m /= m2;
	return m;
}


bool EMAN::operator==(const Matrix3f & m1, const Matrix3f & m2)
{
	const int n = 3;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (m1[i][j] != m2[i][j]) {
				return false;
			}
		}
	}
	return true;
}


bool EMAN::operator!=(const Matrix3f & m1, const Matrix3f & m2)
{
	return !(m1 == m2);
}

Vec3 < float >EMAN::operator*(const Vec3 < int >&v, const Matrix3f & m1)
{
	Vec3 < float >result;
	for (int i = 0; i < 3; i++) {
		result[i] = 0;
		for (int j = 0; j < 3; j++) {
			result[i] += (float)m1[i][j] * v[j];
		}
	}
	return result;
}

Vec3 < float >EMAN::operator*(const Matrix3f & m1, const Vec3 < int >&v)
{
	Vec3 < float >result;
	for (int i = 0; i < 3; i++) {
		result[i] = 0;
		for (int j = 0; j < 3; j++) {
			result[i] += (float)m1[i][j] * v[j];
		}
	}
	return result;
}

Vec3 < float >EMAN::operator*(const Vec3 < float >&v, const Matrix3f & m1)
{
	Vec3 < float >result;
	for (int i = 0; i < 3; i++) {
		result[i] = 0;
		for (int j = 0; j < 3; j++) {
			result[i] += (float) m1[i][j] * v[j];
		}
	}
	return result;
}

Vec3 < float >EMAN::operator*(const Matrix3f & m1, const Vec3 < float >&v)
{
	Vec3 < float >result;
	for (int i = 0; i < 3; i++) {
		result[i] = 0;
		for (int j = 0; j < 3; j++) {
			result[i] += (float) m1[i][j] * v[j];
		}
	}
	return result;
}


/////////////////////////

Matrix4f::Matrix4f()
{
	matrix = gsl_matrix_alloc(4, 4);
	gsl_matrix_set_identity(matrix);
}


Matrix4f::Matrix4f(const vector < float >&m)
{
	const int n = 4;
	matrix = gsl_matrix_alloc(n, n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			gsl_matrix_set(matrix, i, j, m[i * n + j]);
		}
	}
}

Matrix4f::Matrix4f(const Matrix3f & m3)
{
	matrix = gsl_matrix_alloc(4, 4);
	gsl_matrix_set_identity(matrix);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			gsl_matrix_set(matrix, i, j, m3[i][j]);
		}
	}
}


Matrix4f::Matrix4f(gsl_matrix * m)
{
	const int n = 4;
	matrix = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(matrix, m);
}

Matrix4f::~Matrix4f()
{
	if (matrix) {
		gsl_matrix_free(matrix);
		matrix = 0;
	}
}

Matrix4f::Matrix4f(const Matrix4f & m)
{
	const int n = 4;
	matrix = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(matrix, m.get_gsl_matrix());
}


Matrix4f & Matrix4f::operator=(const Matrix4f & m)
{
	if (this != &m) {
		if (matrix) {
			gsl_matrix_free(matrix);
		}
		const int n = 4;
		matrix = gsl_matrix_alloc(n, n);
		gsl_matrix_memcpy(matrix, m.get_gsl_matrix());
	}
	return *this;
}

Matrix4f & Matrix4f::mult_right(const Matrix4f & m)
{
	(*this) = (*this) * m;
	return *this;
}

Matrix4f & Matrix4f::mult_left(const Matrix4f & m)
{
	Matrix4f mcopy = m;
	mcopy = mcopy * (*this);
	(*this) = mcopy;
	return *this;
}

gsl_matrix *Matrix4f::get_gsl_matrix() const
{
	return matrix;
}

void Matrix4f::make_identity()
{
	gsl_matrix_set_identity(matrix);
}


void Matrix4f::set_value(const vector < float >&m)
{
	const int n = 4;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			gsl_matrix_set(matrix, i, j, m[i * n + j]);
		}
	}
}


vector < float >Matrix4f::get_as_list() const
{
	const int n = 4;
	vector < float >m(n * n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			m[i * n + j] = gsl_matrix_get(matrix, i, j);
		}
	}
	return m;
}

Matrix3f Matrix4f::get_matrix3() const
{
	const int n = 3;
	Matrix3f m;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			m[i][j] = gsl_matrix_get(matrix, i, j);
		}
	}
	return m;
}


Matrix4f & Matrix4f::inverse()
{
	const int n = 4;
	gsl_matrix *inverse = gsl_matrix_alloc(n, n);
	gsl_permutation *p = gsl_permutation_alloc(n);

	int s = 0;
	gsl_linalg_LU_decomp(matrix, p, &s);
	gsl_linalg_LU_invert(matrix, p, inverse);
	gsl_matrix_memcpy(matrix, inverse);
	gsl_matrix_free(inverse);
	inverse = 0;

	return (*this);
}

Matrix4f & Matrix4f::transpose()
{
	gsl_matrix_transpose(matrix);
	return (*this);
}

Matrix4f Matrix4f::create_inverse() const
{
	const int n = 4;
	gsl_matrix *matrix_copy = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(matrix_copy, matrix);

	gsl_matrix *inverse = gsl_matrix_alloc(n, n);
	gsl_permutation *p = gsl_permutation_alloc(n);

	int s = 0;
	gsl_linalg_LU_decomp(matrix_copy, p, &s);
	gsl_linalg_LU_invert(matrix_copy, p, inverse);

	gsl_matrix_free(matrix_copy);
	matrix_copy = 0;

	return Matrix4f(inverse);
}


Matrix4f & Matrix4f::operator+=(float f)
{
	gsl_matrix_add_constant(matrix, f);
	return *this;
}

Matrix4f & Matrix4f::operator-=(float f)
{
	gsl_matrix_add_constant(matrix, -f);
	return *this;
}

Matrix4f & Matrix4f::operator*=(float scale)
{
	gsl_matrix_scale(matrix, scale);
	return *this;
}

Matrix4f & Matrix4f::operator/=(float scale)
{
	if (scale != 0) {
		gsl_matrix_scale(matrix, scale);
	}
	return *this;
}


Matrix4f & Matrix4f::operator+=(const Matrix4f & m)
{
	gsl_matrix_add(matrix, m.get_gsl_matrix());
	return *this;
}

Matrix4f & Matrix4f::operator-=(const Matrix4f & m)
{
	gsl_matrix_sub(matrix, m.get_gsl_matrix());
	return *this;
}

Matrix4f & Matrix4f::operator*=(const Matrix4f & m)
{
	int n=4;
	gsl_matrix* temp = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(temp, matrix);
	gsl_matrix_set_zero(matrix);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
				   temp, m.get_gsl_matrix(), 0.0, matrix);
	gsl_matrix_free(temp);
	temp = 0;
	
	return *this;
}

Matrix4f & Matrix4f::operator/=(const Matrix4f & m)
{
	(*this) *= m.create_inverse();
	return *this;
}

double *Matrix4f::operator[] (int i)
{
	return gsl_matrix_ptr(matrix, i, 0);
}

const double *Matrix4f::operator[] (int i) const
{
	return gsl_matrix_const_ptr(matrix, i, 0);
}



Matrix4f EMAN::operator+(float f, const Matrix4f & m2)
{
	Matrix4f m = m2;
	m += f;
	return m;
}

Matrix4f EMAN::operator-(float f, const Matrix4f & m2)
{
	Matrix4f m = m2;
	m -= f;
	return m;
}

Matrix4f EMAN::operator*(float scale, const Matrix4f & m2)
{
	Matrix4f m = m2;
	m *= scale;
	return m;
}

Matrix4f EMAN::operator/(float scale, const Matrix4f & m2)
{
	Matrix4f m = m2;
	m /= scale;
	return m;
}


Matrix4f EMAN::operator+(const Matrix4f & m1, float f)
{
	Matrix4f m = m1;
	m += f;
	return m;
}

Matrix4f EMAN::operator-(const Matrix4f & m1, float f)
{
	Matrix4f m = m1;
	m -= f;
	return m;
}

Matrix4f EMAN::operator*(const Matrix4f & m1, float scale)
{
	Matrix4f m = m1;
	m *= scale;
	return m;
}

Matrix4f EMAN::operator/(const Matrix4f & m1, float scale)
{
	Matrix4f m = m1;
	m /= scale;
	return m;
}


Matrix4f EMAN::operator+(const Matrix4f & m1, const Matrix4f & m2)
{
	Matrix4f m = m1;
	m += m2;
	return m;
}

Matrix4f EMAN::operator-(const Matrix4f & m1, const Matrix4f & m2)
{
	Matrix4f m = m1;
	m -= m2;
	return m;
}

Matrix4f EMAN::operator*(const Matrix4f & m1, const Matrix4f & m2)
{
	Matrix4f m = m1;
	m *= m2;
	return m;
}

Matrix4f EMAN::operator/(const Matrix4f & m1, const Matrix4f & m2)
{
	Matrix4f m = m1;
	m /= m2;
	return m;
}


bool EMAN::operator==(const Matrix4f & m1, const Matrix4f & m2)
{
	const int n = 4;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (m1[i][j] != m2[i][j]) {
				return false;
			}
		}
	}
	return true;
}


bool EMAN::operator!=(const Matrix4f & m1, const Matrix4f & m2)
{
	return !(m1 == m2);
}



/////////////////////////


Quaternion::Quaternion()
:	e0(0), e1(0), e2(0), e3(0)
{
}

Quaternion::Quaternion(float radians, const Vec3 < float >&axis)
{
	Vec3 < float >q = axis;
	//normalize();

	q *= sin(radians / 2.0f);
	e0 = cos(radians / 2.0f);

	e1 = q[0];
	e2 = q[1];
	e3 = q[2];
}

Quaternion::Quaternion(const Vec3 < float >&axis, float radians)
{
	Vec3 < float >q = axis;
	//normalize();

	q *= sin(radians / 2.0f);
	e0 = cos(radians / 2.0f);

	e1 = q[0];
	e2 = q[1];
	e3 = q[2];
}

Quaternion::Quaternion(float ee0, float ee1, float ee2, float ee3)
:e0(ee0), e1(ee1), e2(ee2), e3(ee3)
{
	//normalize();
}

Quaternion::Quaternion(const Matrix3f & m)
{
	int i = 0;

	if (m[0][0] > m[1][1]) {
		if (m[0][0] > m[2][2]) {
			i = 0;
		}
		else {
			i = 2;
		}
	}
	else {
		if (m[1][1] > m[2][2]) {
			i = 1;
		}
		else {
			i = 2;
		}
	}

	if (m[0][0] + m[1][1] + m[2][2] > m[i][i]) {
		e0 = sqrt(m[0][0] + m[1][1] + m[2][2] + 1) / 2.0f;
		e1 = (m[1][2] - m[2][1]) / (4 * e0);
		e2 = (m[2][0] - m[0][2]) / (4 * e0);
		e3 = (m[0][1] - m[1][0]) / (4 * e0);
	}
	else {
		float quat[3];
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;

		quat[i] = sqrt(m[i][i] - m[j][j] - m[k][k] + 1) / 2.0f;
		quat[j] = (m[i][j] + m[j][i]) / (4 * quat[i]);
		quat[k] = (m[i][k] + m[k][i]) / (4 * quat[i]);

		e0 = (m[j][k] - m[k][j]) / (4 * quat[i]);
		e1 = quat[0];
		e2 = quat[1];
		e3 = quat[2];
	}

	//normalize();
}


Quaternion::Quaternion(const Matrix4f & m)
{
	int i = 0;

	if (m[0][0] > m[1][1]) {
		if (m[0][0] > m[2][2]) {
			i = 0;
		}
		else {
			i = 2;
		}
	}
	else {
		if (m[1][1] > m[2][2]) {
			i = 1;
		}
		else {
			i = 2;
		}
	}

	if (m[0][0] + m[1][1] + m[2][2] > m[i][i]) {
		e0 = sqrt(m[0][0] + m[1][1] + m[2][2] + m[3][3]) / 2.0f;
		e1 = (m[1][2] - m[2][1]) / (4 * e0);
		e2 = (m[2][0] - m[0][2]) / (4 * e0);
		e3 = (m[0][1] - m[1][0]) / (4 * e0);
	}
	else {
		float quat[3];
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;

		quat[i] = sqrt(m[i][i] - m[j][j] - m[k][k] + m[3][3]) / 2.0f;
		quat[j] = (m[i][j] + m[j][i]) / (4 * quat[i]);
		quat[k] = (m[i][k] + m[k][i]) / (4 * quat[i]);

		e0 = (m[j][k] - m[k][j]) / (4 * quat[i]);
		e1 = quat[0];
		e2 = quat[1];
		e3 = quat[2];
	}

	//normalize();
}


void Quaternion::normalize()
{
	float dist = 1.0f / sqrt(norm());
	e0 *= dist;
	e1 *= dist;
	e2 *= dist;
	e3 *= dist;
}

Quaternion & Quaternion::inverse()
{
	float f = 1.0f / norm();
	e0 *= f;
	e1 *= -f;
	e2 *= -f;
	e3 *= -f;
	return (*this);
}

Quaternion Quaternion::create_inverse() const
{
	Quaternion q = *this;
	return q.inverse();
}


Vec3 < float >Quaternion::rotate(const Vec3 < float >&v) const
{
	Vec3 < float >i(e1, e2, e3);
	Vec3 < float >v1 = i.cross(v) * (2 * e0);
	Vec3 < float >v2 = v1.cross(i) * (float) 2;
	Vec3 < float >rotated = v + v1 - v2;

	return rotated;
}


float Quaternion::to_angle() const
{
	Vec3 < float >q(e1, e2, e3);
	float len = q.length();
	float radians = 0;

	if (len > 0.00001f) {
		radians = 2.0f * acos(e0);
	}
	else {
		radians = 0;
	}
	return radians;
}

Vec3 < float >Quaternion::to_axis() const
{
	Vec3 < float >q(e1, e2, e3);
	float len = q.length();
	Vec3 < float >axis;

	if (len > 0.00001f) {
		axis = q * ((float) (1.0f / len));
	}
	else {
		axis.set_value(0.0f, 0.0f, 1.0f);
	}
	return axis;
}


Matrix3f Quaternion::to_matrix3() const
{
	vector < float >m(9);

	m[0] = e0 * e0 + e1 * e1 - e2 * e2 - e3 * e3;
	m[1] = 2.0f * (e1 * e2 + e0 * e3);
	m[2] = 2.0f * (e1 * e3 - e0 * e2);

	m[3] = 2.0f * (e1 * e2 - e0 * e3);
	m[4] = e0 * e0 + e1 * e1 + e2 * e2 - e3 * e3;
	m[5] = 2.0f * (e2 * e3 + e0 * e1);

	m[6] = 2.0f * (e1 * e3 + e0 * e2);
	m[7] = 2.0f * (e2 * e3 - e0 * e1);
	m[8] = e0 * e0 - e1 * e1 - e2 * e2 + e3 * e3;

	return Matrix3f(m);
}

Matrix4f Quaternion::to_matrix4() const
{
	vector < float >m(16);

	m[0] = e0 * e0 + e1 * e1 - e2 * e2 - e3 * e3;
	m[1] = 2.0f * (e1 * e2 + e0 * e3);
	m[2] = 2.0f * (e1 * e3 - e0 * e2);
	m[3] = 0;

	m[4] = 2.0f * (e1 * e2 - e0 * e3);
	m[5] = e0 * e0 + e1 * e1 + e2 * e2 - e3 * e3;
	m[6] = 2.0f * (e2 * e3 + e0 * e1);
	m[7] = 0;

	m[8] = 2.0f * (e1 * e3 + e0 * e2);
	m[9] = 2.0f * (e2 * e3 - e0 * e1);
	m[10] = e0 * e0 - e1 * e1 - e2 * e2 + e3 * e3;
	m[11] = 0;

	m[12] = 0;
	m[13] = 0;
	m[14] = 0;
	m[15] = e0 * e0 + e1 * e1 + e2 * e2 + e3 * e3;

	return Matrix4f(m);
}



float Quaternion::real() const
{
	return e0;
}

Vec3 < float >Quaternion::unreal() const
{
	return Vec3 < float >(e1, e2, e3);
}

vector < float >Quaternion::get_as_list() const
{
	vector < float >v(4);
	v[0] = e0;
	v[1] = e1;
	v[2] = e2;
	v[3] = e3;

	return v;
}

Quaternion & Quaternion::operator+=(const Quaternion & q)
{
	e0 += q.e0;
	e1 += q.e1;
	e2 += q.e2;
	e3 += q.e3;
	return *this;
}

Quaternion & Quaternion::operator-=(const Quaternion & q)
{
	e0 -= q.e0;
	e1 -= q.e1;
	e2 -= q.e2;
	e3 -= q.e3;
	return *this;
}

Quaternion & Quaternion::operator*=(const Quaternion & q)
{
	float a = e0 * q.e0 - e1 * q.e1 - e2 * q.e2 - e3 * q.e3;
	float b = e0 * q.e1 + e1 * q.e0 + e2 * q.e3 - e3 * q.e2;
	float c = e0 * q.e2 - e1 * q.e3 + e2 * q.e0 + e3 * q.e1;
	float d = e0 * q.e3 + e1 * q.e2 - e2 * q.e1 + e3 * q.e0;

	e0 = a;
	e1 = b;
	e2 = c;
	e3 = d;


	// normalize();

	return (*this);
}

Quaternion & Quaternion::operator*=(float s)
{
	e0 *= s;
	e1 *= s;
	e2 *= s;
	e3 *= s;
	return (*this);
}

Quaternion & Quaternion::operator/=(const Quaternion & q)
{
	float qn = q.norm();

	float a = (+e0 * q.e0 + e1 * q.e1 + e2 * q.e2 + e3 * q.e3) / qn;
	float b = (-e0 * q.e1 + e1 * q.e0 - e2 * q.e3 + e3 * q.e2) / qn;
	float c = (-e0 * q.e2 + e1 * q.e3 + e2 * q.e0 - e3 * q.e1) / qn;
	float d = (-e0 * q.e3 - e1 * q.e2 + e2 * q.e1 + e3 * q.e0) / qn;

	e0 = a;
	e1 = b;
	e2 = c;
	e3 = d;

	return (*this);
}

Quaternion & Quaternion::operator/=(float s)
{
	if (s != 0) {
		e0 /= s;
		e1 /= s;
		e2 /= s;
		e3 /= s;
	}

	return (*this);
}


Quaternion EMAN::operator+(const Quaternion & q1, const Quaternion & q2)
{
	Quaternion q = q1;
	q += q2;
	return q;
}

Quaternion EMAN::operator-(const Quaternion & q1, const Quaternion & q2)
{
	Quaternion q = q1;
	q -= q2;
	return q;
}


Quaternion EMAN::operator*(const Quaternion & q1, const Quaternion & q2)
{
	Quaternion q = q1;
	q *= q2;
	return q;
}

Quaternion EMAN::operator*(const Quaternion & q, float s)
{
	Quaternion q1 = q;
	q1 *= s;
	return q1;
}

Quaternion EMAN::operator*(float s, const Quaternion & q)
{
	Quaternion q1 = q;
	q1 *= s;
	return q1;
}

Quaternion EMAN::operator/(const Quaternion & q1, const Quaternion & q2)
{
	Quaternion q = q1;
	q /= q2;
	return q;
}


bool EMAN::operator==(const Quaternion & q1, const Quaternion & q2)
{
	bool result = true;
	vector < float >v1 = q1.get_as_list();
	vector < float >v2 = q2.get_as_list();

	for (size_t i = 0; i < v1.size(); i++) {
		if (v1[i] != v2[i]) {
			result = false;
			break;
		}
	}

	return result;
}

bool EMAN::operator!=(const Quaternion & q1, const Quaternion & q2)
{
	return (!(q1 == q2));
}


Quaternion Quaternion::interpolate(const Quaternion & from, const Quaternion & to, float t)
{
	const double epsilon = 0.00001f;
	double cosom = from.e1 * to.e1 + from.e2 * to.e2 + from.e3 * to.e3 + from.e0 * to.e0;

	Quaternion q;
	if (cosom < 0.0f) {
		cosom = -cosom;
		q = q - to;
	}
	else {
		q = to;
	}

	double scale0 = 1 - t;
	double scale1 = t;

	if ((1 - cosom) > epsilon) {
		double omega = acos(cosom);
		double sinom = sin(omega);
		scale0 = sin((1 - t) * omega) / sinom;
		scale1 = sin(t * omega) / sinom;
	}

	return (scale0 * from + scale1 * q);
}


/////////////////////////


const float Rotation::ERR_LIMIT = 0.000001;


Rotation::Rotation()
:	type(UNKNOWN), symname(""), nsym(0), cur_sym(0)
{
	matrix.make_identity();
}


Rotation::Rotation(float a1, float a2, float a3, Type t)
:type(t), symname(""), nsym(0), cur_sym(0)
{
	set_angle(a1, a2, a3, t);
}


Rotation::Rotation(float e1, float e2, float e3, float e0, Type t)
:type(t), symname(""), nsym(0), cur_sym(0)
{
	set_angle(e1, e2, e3, e0, t);
}

Rotation::Rotation(const Quaternion & q)
:type(QUATERNION), symname(""), nsym(0), cur_sym(0)
{
	set_angle(q);
}

Rotation::Rotation(const Matrix3f & m)
:matrix(m), type(MATRIX), symname(""), nsym(0), cur_sym(0)
{
}

Rotation::~Rotation()
{
}


Rotation & Rotation::inverse()
{
	matrix.inverse();
	return (*this);
}

Rotation Rotation::create_inverse() const
{
	return Rotation(matrix.create_inverse());
}


float Rotation::diff(const Rotation & r)
{
	float denom = sqrt((float) 3);
	Matrix3f matrix2 = r.get_matrix3();

	float a1 = (matrix[0][0] + matrix[0][1] + matrix[0][2]) / denom;
	float a2 = (matrix[1][0] + matrix[1][1] + matrix[1][2]) / denom;
	float a3 = (matrix[2][0] + matrix[2][1] + matrix[2][2]) / denom;

	float b1 = (matrix2[0][0] + matrix2[0][1] + matrix2[0][2]) / denom;
	float b2 = (matrix2[1][0] + matrix2[1][1] + matrix2[1][2]) / denom;
	float b3 = (matrix2[2][0] + matrix2[2][1] + matrix2[2][2]) / denom;

	float p = (a1 * b1 + a2 * b2 + a3 * b3);
	if (p > 1) {
		p = 1;
	}
	else if (p < -1) {
		p = -1;
	}
	float result = acos(p);
	return result;
}

Rotation & Rotation::rotate_from_left(const Rotation & left)
{
	Matrix3f m = left.get_matrix3();
	m = m * matrix;
	matrix = m;
	return (*this);
}

Rotation::SymType Rotation::get_sym_type(string name)
{
	SymType t = UNKNOWN_SYM;

	if (name[0] == 'c') {
		t = CSYM;
	}
	else if (name[0] == 'd') {
		t = DSYM;
	}
	else if (name == "icos") {
		t = ICOS_SYM;
	}
	else if (name == "oct") {
		t = OCT_SYM;
	}
	else if (name == "i") {
		t = ISYM;
	}
	return t;
}


void Rotation::set_sym(string name)
{
	if (name == "") {
		symname = "i";
		return;
	}

	symname = name;

	for (size_t i = 0; i < name.size(); i++) {
		if (isalpha(name[i])) {
			symname[i] = tolower(name[i]);
		}
	}

	SymType type = get_sym_type(symname);
	switch (type) {
	case CSYM:
		nsym = atoi(symname.c_str() + 1);
		break;
	case DSYM:
		nsym = atoi(symname.c_str() + 1) * 2;
		break;
	case ICOS_SYM:
		nsym = 60;
		break;
	case OCT_SYM:
		nsym = 24;
		break;
	case ISYM:
		nsym = 1;
		break;
	case UNKNOWN_SYM:
	default:
		nsym = 0;
	}
}


int Rotation::get_max_nsym() const
{
	return nsym;
}

Rotation Rotation::get_sym(int n)
{
	Rotation invalid(-0.1, -0.1, -0.1, Rotation::EMAN);
	Rotation ret;

	static float ICOS[180] = {
		0, 0, 0, 0, 0, 288, 0, 0, 216, 0, 0, 144, 0, 0, 72,
		0, 63.4349, 36, 0, 63.4349, 324, 0, 63.4349, 252, 0, 63.4349, 180, 0, 63.4349, 108,
		72, 63.4349, 36, 72, 63.4349, 324, 72, 63.4349, 252, 72, 63.4349, 180, 72, 63.4349, 108,
		144, 63.4349, 36, 144, 63.4349, 324, 144, 63.4349, 252, 144, 63.4349, 180, 144, 63.4349,
		108,
		216, 63.4349, 36, 216, 63.4349, 324, 216, 63.4349, 252, 216, 63.4349, 180, 216, 63.4349,
		108,
		288, 63.4349, 36, 288, 63.4349, 324, 288, 63.4349, 252, 288, 63.4349, 180, 288, 63.4349,
		108,
		36, 116.5651, 0, 36, 116.5651, 288, 36, 116.5651, 216, 36, 116.5651, 144, 36, 116.5651, 72,
		108, 116.5651, 0, 108, 116.5651, 288, 108, 116.5651, 216, 108, 116.5651, 144, 108, 116.5651,
		72,
		180, 116.5651, 0, 180, 116.5651, 288, 180, 116.5651, 216, 180, 116.5651, 144, 180, 116.5651,
		72,
		252, 116.5651, 0, 252, 116.5651, 288, 252, 116.5651, 216, 252, 116.5651, 144, 252, 116.5651,
		72,
		324, 116.5651, 0, 324, 116.5651, 288, 324, 116.5651, 216, 324, 116.5651, 144, 324, 116.5651,
		72,
		0, 180, 0, 0, 180, 288, 0, 180, 216, 0, 180, 144, 0, 180, 72
	};

	static float OCT[72] = {
		0, 0, 0, 90, 0, 0, 180, 0, 0, 270, 0, 0, 90, 90, 0, 90, 270, 0,
		0, 0, 90, 90, 0, 90, 180, 0, 90, 270, 0, 90, 90, 90, 90, 90, 270, 90,
		0, 0, 180, 90, 0, 180, 180, 0, 180, 270, 0, 180, 90, 90, 180, 90, 270, 180,
		0, 0, 270, 90, 0, 270, 180, 0, 270, 270, 0, 270, 90, 90, 270, 90, 270, 270
	};

	if (n >= nsym) {
		return invalid;
	}

	SymType type = get_sym_type(symname);

	switch (type) {
	case CSYM:
		ret.set_angle(0, n * 2.0f * M_PI / nsym, 0, Rotation::EMAN);
		break;
	case DSYM:
		if (n >= nsym / 2) {
			ret.set_angle(M_PI, (n - nsym / 2) * 2.0f * M_PI / (nsym / 2), 0, Rotation::EMAN);
		}
		else {
			ret.set_angle(0, n * 2.0f * M_PI / (nsym / 2), 0, Rotation::EMAN);
		}
		break;
	case ICOS_SYM:
		ret.set_angle(ICOS[n * 3 + 1] * M_PI / 180.0f, ICOS[n * 3 + 2] * M_PI / 180.0f - M_PI / 2,
					  ICOS[n * 3] * M_PI / 180.0f + M_PI / 2., Rotation::EMAN);
		break;
	case OCT_SYM:
		ret.set_angle(OCT[n * 3] * M_PI / 180.0f, OCT[n * 3 + 1] * M_PI / 180.0f,
					  OCT[n * 3 + 2] * M_PI / 180.0f, Rotation::EMAN);
		break;
	case ISYM:
		ret.set_angle(0, 0, 0, Rotation::EMAN);
	default:
		LOGERR("I don't know this symmetry: %s", symname.c_str());
	}

	ret = (*this) * ret;

	return ret;
}

void Rotation::rectify()
{

}

/* x-convention: z-axis, x-axis, z-axis.
   
in EMAN convention, "x-convention" is used, but a1 is second angle (alt),
a2 is the first angle (az), a3 is the 3rd (phi)
MRC convention is 'y-convention'.

The Euler angles in SPIDER & Web are defined as three successive
rotations in a right hand coordinate system. First, the object is
rotated COUNTERCLOCKWISE around the Z-axis (angle 'phi') and then it
is rotated CLOCKWISE around the Y-axis (angle 'theta') and finally, it
is roated COUNTERCLOCKWISE around the Z-axis (angle 'psi'). All
rotations are done around axes of the SPACE coordinate system.
*/

void Rotation::set_angle(float a1, float a2, float a3, Type t)
{
	type = t;
	if (type == IMAGIC) {
		// nothing
	}
	else if (type == SPIDER) {
		a2 = a2 - 1.5 * M_PI;
		a3 = a3 - 0.5 * M_PI;
	}
	else if (type == MRC) {
		a2 = fmod(-a2 + 2.5 * M_PI, 2 * M_PI);
		a3 = fmod(a3 + 0.5 * M_PI, 2.0f * M_PI);
	}
	else if (type == EMAN) {
		// nothing
	}
	else {
		LOGERR("unknown Euler angle convention: %d\n", type);
		return;
	}
	// for eman: a1: alt, a2: az, a3: phi.
#if 0
	// this version is from EMAN1 Euler.C; seems not working
	matrix[0][0] = cos(a3) * cos(a2) - cos(a1) * sin(a2) * sin(a3);
	matrix[0][1] = cos(a3) * sin(a2) + cos(a1) * cos(a2) * sin(a3);
	matrix[0][2] = sin(a3) * sin(a1);

	matrix[1][0] = -sin(a3) * cos(a2) - cos(a1) * sin(a2) * cos(a3);
	matrix[1][1] = -sin(a3) * sin(a2) + cos(a1) * cos(a2) * cos(a3);
	matrix[1][2] = cos(a3) * sin(a1);

	matrix[2][0] = sin(a1) * sin(a2);
	matrix[2][1] = -sin(a1) * cos(a2);
	matrix[2][2] = cos(a1);
#endif

	// the following version is from EMAN1 EMDataB.C rotateAndTranslate().

	matrix[0][0] = cos(a3)*cos(a2) - cos(a1)*sin(a2)*sin(a3);
	matrix[0][1] = -(sin(a3)*cos(a2) + cos(a1)*sin(a2)*cos(a3));
	matrix[0][2] = sin(a1)*sin(a2);
	matrix[1][0] = cos(a3)*sin(a2) + cos(a1)*cos(a2)*sin(a3);
	matrix[1][1] = -sin(a3)*sin(a2) + cos(a1)*cos(a2)*cos(a3);
	matrix[1][2] = -sin(a1)*cos(a2);
	matrix[2][0] = sin(a1)*sin(a3);
	matrix[2][1] = sin(a1)*cos(a3);
	matrix[2][2] = cos(a1);

	
	rectify();
}

void Rotation::set_angle(float e1, float e2, float e3, float e0, Type t)
{
	type = t;
	if (type != SPIN && type != SGIROT && type != QUATERNION) {
		LOGERR("unknown Euler angle convention: %d\n", type);
		return;
	}

	Quaternion *q = 0;
	if (type == QUATERNION) {
		q = new Quaternion(e0, e1, e2, e3);
	}
	else {
		q = new Quaternion(e0, Vec3 < float >(e1, e2, e3));
	}

	set_angle(*q);

	if (type == SGIROT) {
		inverse();
	}

	delete q;
	q = 0;
}

void Rotation::set_angle(const Matrix3f & m)
{
	matrix = m;
}

void Rotation::set_angle(const Quaternion & q)
{
	matrix = q.to_matrix3();
}


bool Rotation::is_valid() const
{
	return true;
}


float Rotation::eman_alt() const
{
	float alt = 0;
	float max = 1 - ERR_LIMIT;

	if (matrix[2][2] > max) {
		alt = 0;
	}
	else if (matrix[2][2] < -max) {
		alt = M_PI;
	}
	else {
		alt = acos(matrix[2][2]);
	}
	return alt;
}


float Rotation::eman_az() const
{
	float az = 0;
	float max = 1 - ERR_LIMIT;
	if (matrix[2][2] > max) {
		az = atan2(matrix[0][1], matrix[1][1]);
	}
	else if (matrix[2][2] < -max) {
		az = atan2(matrix[0][1], -matrix[1][1]);
	}
	else {
		az = atan2(matrix[2][0], -matrix[2][1]);
	}
	return az;
}

float Rotation::eman_phi() const
{
	float phi = 0;

	if (fabs(matrix[2][2]) > (1 - ERR_LIMIT)) {
		phi = 0;
	}
	else {
		phi = atan2(matrix[0][2], matrix[1][2]);
	}

	return phi;
}


float Rotation::mrc_theta() const
{
	return eman_alt();
}

float Rotation::mrc_phi() const
{
	float az = eman_az();
	return fmod(-az + 2.5 * M_PI, 2 * M_PI);
}

float Rotation::mrc_omega() const
{
	float phi = eman_phi();
	return fmod(phi + 1.5 * M_PI, 2 * M_PI);
}


float Rotation::imagic_alpha() const
{
	return eman_phi();
}

float Rotation::imagic_beta() const
{
	return eman_alt();
}

float Rotation::imagic_gamma() const
{
	return eman_az();
}


float Rotation::spider_phi() const
{
	return fmod((eman_az() + 1.5 * M_PI), 2.0 * M_PI);
}

float Rotation::spider_theta() const
{
	return (eman_alt());
}

float Rotation::spider_gamma() const
{
	return fmod((eman_phi() + 0.5 * M_PI), 2.0 * M_PI);
}

vector < float >Rotation::get_spin_axis() const
{
	Quaternion quat(matrix);
	float q = quat.to_angle();
	Vec3 < float >v = quat.to_axis();

	vector < float >result(4);
	result[0] = q;
	result[1] = v[0];
	result[2] = v[1];
	result[3] = v[2];

	return result;
}

void Rotation::get_spin_axis(float *q, float *n1, float *n2, float *n3) const
{
	Quaternion quat(matrix);
	*q = quat.to_angle();
	Vec3 < float >v = quat.to_axis();
	*n1 = v[0];
	*n2 = v[1];
	*n3 = v[2];
}

vector < float >Rotation::get_sgi() const
{
	return get_spin_axis();
}

void Rotation::get_sgi(float *q, float *n1, float *n2, float *n3) const
{
	get_spin_axis(q, n1, n2, n3);
}


Quaternion Rotation::get_quaternion() const
{
	return Quaternion(matrix);
}


Matrix3f Rotation::get_matrix3() const
{
	return matrix;
}

Matrix4f Rotation::get_matrix4() const
{
	return Matrix4f(matrix);
}

Rotation & Rotation::operator*=(const Rotation & e)
{
	matrix *= e.get_matrix3();
	return (*this);
}

Rotation & Rotation::operator/=(const Rotation & e)
{
	Matrix3f m2 = e.get_matrix3();
	m2.inverse();
	matrix *= m2;
	return (*this);
}


Rotation EMAN::operator*(const Rotation & e1, const Rotation & e2)
{
	Rotation r = e1;
	r *= e2;
	return r;
}

Rotation EMAN::operator/(const Rotation & e1, const Rotation & e2)
{
	Rotation r = e1;
	r /= e2;
	return r;
}


bool EMAN::operator==(const Rotation & e1, const Rotation & e2)
{
	return (e1.get_matrix3() == e2.get_matrix3());
}

bool EMAN::operator!=(const Rotation & e1, const Rotation & e2)
{
	return (!(e1 == e2));
}


//////////////////////


Transform & Transform::set_transform_instance(const Vec3 < float >&translation,
											  const Rotation & rotation,
											  const Vec3 < float >&scale_factor,
											  const Rotation & scale_orientation,
											  const Vec3 < float >&center)
{
	Rotation so = scale_orientation;

	matrix.make_identity();
	translate((-1.0f) * center);

	if (so != Rotation(1, 0, 0, 0, Rotation::QUATERNION)) {
		Rotation inverse_so = so.create_inverse();
		rotate(inverse_so);
	}
	scale(scale_factor);
	rotate(so);
	rotate(rotation);
	translate(center);
	translate(translation);

	return (*this);
}



Vec3 < float >Transform::transform(const Vec3 < float >&x)
{
	Vec3 < float >x1 = x - pre_trans;
	Matrix3f r = matrix.get_matrix3();
	x1 = r * x1;
	x1 += get_post_translate();

	return x1;
}


Vec3 < float >Transform::inverse_transform(const Vec3 < float >&v)
{
	Transform t = create_inverse();
	return t.transform(v);
}

float Transform::get_scale(int i) const
{
	float s = matrix[i][0] *  matrix[i][0] + matrix[i][1]*matrix[i][1]
		+ matrix[i][2] * matrix[i][2];
	return s;
}

Vec3 < float > Transform::get_scale() const
{
	float s1 = matrix[0][0] *  matrix[0][0] + matrix[0][1]*matrix[0][1]
		+ matrix[0][2] * matrix[0][2];
	float s2 = matrix[1][0] *  matrix[1][0] + matrix[1][1]*matrix[1][1]
		+ matrix[1][2] * matrix[1][2];

	float s3 = matrix[2][0] *  matrix[2][0] + matrix[2][1]*matrix[2][1]
		+ matrix[2][2] * matrix[2][2];
	
	return Vec3 < float >(s1, s2, s3);
}


int Transform::get_type() const
{
	return TRANSFORM;
}



Transform EMAN::operator+(const Transform & t1, const Transform & t2)
{
	Transform t = t1;
	t += t2;
	return t;
}

Transform EMAN::operator-(const Transform & t1, const Transform & t2)
{
	Transform t = t1;
	t -= t2;
	return t;
}


Transform EMAN::operator*(const Transform & t1, const Transform & t2)
{
	Transform t = t1;
	t *= t2;
	return t;
}

Transform EMAN::operator*(const Transform & t, float s)
{
	Transform new_t = t;
	new_t *= s;
	return new_t;
}

Transform EMAN::operator*(float s, const Transform & t)
{
	Transform new_t = t;
	new_t *= s;
	return new_t;
}


Transform EMAN::operator/(const Transform & t1, const Transform & t2)
{
	Transform t = t1;
	t /= t2;
	return t;
}

Transform EMAN::operator/(const Transform & t, float s)
{
	Transform new_t = t;
	new_t /= s;
	return new_t;
}

Transform EMAN::operator/(float s, const Transform & t)
{
	Transform new_t = t;
	new_t /= s;
	return new_t;
}

Rotation Rotation::interpolate(const Rotation & from, const Rotation & to, float percent)
{
	Quaternion q1 = from.get_quaternion();
	Quaternion q2 = to.get_quaternion();
	Quaternion q = Quaternion::interpolate(q1, q2, percent);
	return Rotation(q);
}

