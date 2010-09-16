/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * */

#include "transform.h"
#include "util.h"
#include "emobject.h"
#include <cctype> // for std::tolower
#include <cstring>  // for memcpy
#include "symmetry.h"
using namespace EMAN;

#ifdef WIN32
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif

#include <algorithm> // for std::transform

#include <gsl_matrix.h>
#include <gsl_blas.h>
#include <gsl_linalg.h>

#include <ostream>
using std::ostream_iterator;

const float Transform3D::ERR_LIMIT = 0.000001f;

const float Transform::ERR_LIMIT = 0.000001f;

vector<string>  Transform::permissable_2d_not_rot;
vector<string>  Transform::permissable_3d_not_rot;
map<string,vector<string> >  Transform::permissable_rot_keys;

Transform Transform::icos_5_to_2() {
	Transform t;
	Dict d;
	d["type"] = "eman";
	d["phi"] = 0.0f;
	d["az"] = 270.0f;
	d["alt"] = 58.282523f; // 5 fold to a 3 fold to a 2 fold
	/** Doctor Phil says:
	 * alt = (acos(cos(pi/5)/sqrt(3)/sin(pi/5)) + acos(2*cos(pi/5)/ sqrt(3) ) )*180/pi
	 * This is the angle between a 5 and a 3 plus the angle between a 3 and a 2
	 */
	
	t.set_rotation(d);
	return t;
}

Transform Transform::tet_3_to_2() {
	Transform t;
	Dict d;
	d["type"] = "eman";
	d["phi"] = 45.0f;
	d["az"] = 0.0f;
	d["alt"] = 54.73561f; // 3 fold to a 2 fold
	/** Doctor Phil says:
	 * AltAngle= acos(-1/3.0)*90/pi;
	 */
	t.set_rotation(d);
	return t;
}


Transform::Transform()
{
	to_identity();
}

Transform::Transform( const Transform& that )
{
	*this = that;
}

Transform& Transform::operator=(const Transform& that ) {

	if (this != &that ) {
		memcpy(matrix,that.matrix,12*sizeof(float));
// 		transform_type = that.transform_type;
	}
	return *this;
}

Transform::Transform(const Dict& d)  {
	to_identity();
	set_params(d);
}


Transform::Transform(const float array[12]) {
	memcpy(matrix,array,12*sizeof(float));
}

Transform::Transform(const vector<float> array)
{
	set_matrix(array);
}

void Transform::set_matrix(const vector<float>& v)
{
	if (v.size() != 12 ) throw InvalidParameterException("The construction array must be of size 12");

	for(int i=0; i<3; ++i) {
		for(int j=0; j<4; ++j) {
			matrix[i][j] = v[i*4+j];
		}
	}
}

void Transform::copy_matrix_into_array(float* const array) const {

	int idx = 0;
	for(int i=0; i<3; ++i) {
		for(int j=0; j<4; ++j) {
			array[idx] = matrix[i][j];
			idx ++;
		}
	}
}

vector<float> Transform::get_matrix() const
{
	vector<float> ret(12);
	for(int i=0; i<3; ++i) {
		for(int j=0; j<4; ++j) {
			ret[i*4+j] = matrix[i][j];
		}
	}
	return ret;
}

void Transform::to_identity()
{
// 	transform_type = UNKNOWN;
	for(int i=0; i<3; ++i) {
		for(int j=0; j<4; ++j) {
			if(i==j) {
				matrix[i][j] = 1;
			}
			else {
				matrix[i][j] = 0;
			}
		}
	}
}

bool Transform::is_identity() const {
	for(int i=0; i<3; ++i) {
		for(int j=0; j<4; ++j) {
			float c = matrix[i][j];
			Util::apply_precision(c,ERR_LIMIT);
			if(i==j) {
				if (c != 1.0) return false;
			}
			else {
				if (c != 0.0) return false;
			}
		}
	}
	return true;
}


void Transform::set_params(const Dict& d) {
	detect_problem_keys(d);

	if (d.has_key_ci("type") ) set_rotation(d);

	if (d.has_key_ci("scale")) {
		float scale = static_cast<float>(d.get_ci("scale"));
		set_scale(scale);
	}

	float dx=0,dy=0,dz=0;

	if (d.has_key_ci("tx")) dx = static_cast<float>(d.get_ci("tx"));
	if (d.has_key_ci("ty")) dy = static_cast<float>(d.get_ci("ty"));
	if (d.has_key_ci("tz")) dz = static_cast<float>(d.get_ci("tz"));

	if ( dx != 0.0 || dy != 0.0 || dz != 0.0 ) {
		set_trans(dx,dy,dz);
	}

	if (d.has_key_ci("mirror")) {
		EMObject e = d.get_ci("mirror");
		if ( (e.get_type() != EMObject::BOOL ) && (e.get_type() != EMObject::INT ) && (e.get_type() != EMObject::UNSIGNEDINT ) )
			throw InvalidParameterException("Error, mirror must be a bool or an int");

		bool mirror = static_cast<bool>(e);
		set_mirror(mirror);
	}
}


void Transform::init_permissable_keys()
{

	permissable_2d_not_rot.push_back("tx");
	permissable_2d_not_rot.push_back("ty");
	permissable_2d_not_rot.push_back("scale");
	permissable_2d_not_rot.push_back("mirror");
	permissable_2d_not_rot.push_back("type");

	permissable_3d_not_rot.push_back("tx");
	permissable_3d_not_rot.push_back("ty");
	permissable_3d_not_rot.push_back("tz");
	permissable_3d_not_rot.push_back("scale");
	permissable_3d_not_rot.push_back("mirror");
	permissable_3d_not_rot.push_back("type");

	vector<string> tmp;
	tmp.push_back("alpha");
	permissable_rot_keys["2d"] = tmp;

	tmp.clear();
	tmp.push_back("alt");
	tmp.push_back("az");
	tmp.push_back("phi");
	permissable_rot_keys["eman"] = tmp;

	tmp.clear();
	tmp.push_back("psi");
	tmp.push_back("theta");
	tmp.push_back("phi");
	permissable_rot_keys["spider"] = tmp;

	tmp.clear();
	tmp.push_back("alpha");
	tmp.push_back("beta");
	tmp.push_back("gamma");
	permissable_rot_keys["imagic"] = tmp;

	tmp.clear();
	tmp.push_back("ztilt");
	tmp.push_back("xtilt");
	tmp.push_back("ytilt");
	permissable_rot_keys["xyz"] = tmp;

	tmp.clear();
	tmp.push_back("phi");
	tmp.push_back("theta");
	tmp.push_back("omega");
	permissable_rot_keys["mrc"] = tmp;

	tmp.clear();
	tmp.push_back("e0");
	tmp.push_back("e1");
	tmp.push_back("e2");
	tmp.push_back("e3");
	permissable_rot_keys["quaternion"] = tmp;

	tmp.clear();
	tmp.push_back("n1");
	tmp.push_back("n2");
	tmp.push_back("n3");
	tmp.push_back("Omega");
	permissable_rot_keys["spin"] = tmp;

	tmp.clear();
	tmp.push_back("n1");
	tmp.push_back("n2");
	tmp.push_back("n3");
	tmp.push_back("q");
	permissable_rot_keys["sgirot"] = tmp;

	tmp.clear();
	tmp.push_back("m11");
	tmp.push_back("m12");
	tmp.push_back("m13");
	tmp.push_back("m21");
	tmp.push_back("m22");
	tmp.push_back("m23");
	tmp.push_back("m31");
	tmp.push_back("m32");
	tmp.push_back("m33");
	permissable_rot_keys["matrix"] = tmp;
}


void Transform::detect_problem_keys(const Dict& d) {
	if (permissable_rot_keys.size() == 0 ) {
		init_permissable_keys();
	}

	vector<string> verification;
	vector<string> problem_keys;
	bool is_2d = false;
	if (d.has_key_ci("type") ) {
		string type = Util::str_to_lower((string)d["type"]);
		bool problem = false;
		if (permissable_rot_keys.find(type) == permissable_rot_keys.end() ) {
			problem_keys.push_back(type);
			problem = true;
		}
		if ( !problem ) {
			vector<string> perm = permissable_rot_keys[type];
			std::copy(perm.begin(),perm.end(),back_inserter(verification));

			if ( type == "2d" ) {
				is_2d = true;
				std::copy(permissable_2d_not_rot.begin(),permissable_2d_not_rot.end(),back_inserter(verification));
			}
		}
	}
	if ( !is_2d ) {
		std::copy(permissable_3d_not_rot.begin(),permissable_3d_not_rot.end(),back_inserter(verification));
	}

	for (Dict::const_iterator it = d.begin(); it != d.end();  ++it) {
		if ( std::find(verification.begin(),verification.end(), it->first) == verification.end() ) {
			problem_keys.push_back(it->first);
		}
	}

	if (problem_keys.size() != 0 ) {
		string error;
		if (problem_keys.size() == 1) {
			error = "Transform Error: The \"" +problem_keys[0]+ "\" key is unsupported";
		} else {
			error = "Transform Error: The ";
			for(vector<string>::const_iterator cit = problem_keys.begin(); cit != problem_keys.end(); ++cit ) {
				if ( cit != problem_keys.begin() ) {
					if (cit == (problem_keys.end() -1) ) error += " and ";
					else error += ", ";
				}
				error += "\"";
				error += *cit;
				error += "\"";
			}
			error += " keys are unsupported";
		}
		throw InvalidParameterException(error);
	}
}

void Transform::set_params_inverse(const Dict& d) {
	detect_problem_keys(d);

	if (d.has_key_ci("type") ) set_rotation(d);

	float dx=0,dy=0,dz=0;
	if (d.has_key_ci("tx")) dx = static_cast<float>(d.get_ci("tx"));
	if (d.has_key_ci("ty")) dy = static_cast<float>(d.get_ci("ty"));
	if (d.has_key_ci("tz")) dz = static_cast<float>(d.get_ci("tz"));

	if ( (dx != 0.0 || dy != 0.0 || dz != 0.0) && d.has_key_ci("type") ) {
		Transform pre_trans;
		pre_trans.set_trans(dx,dy,dz);

		Transform tmp;
		tmp.set_rotation(d);

		if (d.has_key_ci("scale")) {
			float scale = static_cast<float>(d.get_ci("scale"));
			tmp.set_scale(scale);
		}

		Transform solution_trans = tmp*pre_trans;

		if (d.has_key_ci("scale")) {
			Transform tmp;
			float scale = static_cast<float>(d.get_ci("scale"));
			tmp.set_scale(scale);
			solution_trans = solution_trans*tmp;
		}

		tmp = Transform();
		tmp.set_rotation(d);
		solution_trans = solution_trans*tmp;
		set_trans(solution_trans.get_trans());
	}

	if (d.has_key_ci("scale")) {
		float scale = static_cast<float>(d.get_ci("scale"));
		set_scale(scale);
	}

	if (d.has_key_ci("mirror")) {
		EMObject e = d.get_ci("mirror");
		if ( (e.get_type() != EMObject::BOOL ) && (e.get_type() != EMObject::INT ) && (e.get_type() != EMObject::UNSIGNEDINT ) )
			throw InvalidParameterException("Error, mirror must be a bool or an int");

		bool mirror = static_cast<bool>(e);
		set_mirror(mirror);
	}
	invert();
}


Dict Transform::get_params(const string& euler_type) const {
	Dict params = get_rotation(euler_type);

	Vec3f v = get_trans();
	params["tx"] = v[0]; params["ty"] = v[1];

	string type = Util::str_to_lower(euler_type);
	if ( type != "2d") params["tz"] = v[2];

	float scale = get_scale();
	params["scale"] = scale;

	bool mirror = get_mirror();
	params["mirror"] = mirror;

	return params;
}



Dict Transform::get_params_inverse(const string& euler_type) const {
	Transform inv(inverse());

	Dict params = inv.get_rotation(euler_type);
	Vec3f v = inv.get_pre_trans();
	params["tx"] = v[0]; params["ty"] = v[1];

	string type = Util::str_to_lower(euler_type);
	if ( type != "2d") params["tz"] = v[2];

	float scale = inv.get_scale();
	params["scale"] = scale;

	bool mirror = inv.get_mirror();
	params["mirror"] = mirror;

	return params;
}

// Dict Transform::get_params_2d() {
// 	Dict params = get_rotation("2d");
//
// 	Vec2f v = get_trans_2d();
// 	params["tx"] = v[0]; params["ty"] = v[1];
//
// 	float scale = get_scale();
// 	params["scale"] = scale;
//
// 	bool mirror = get_mirror();
// 	params["mirror"] = mirror;
//
// 	return params;
// }


void Transform::set_rotation(const Dict& rotation)
{
	detect_problem_keys(rotation);
	string euler_type;

	if (!rotation.has_key_ci("type") ){
			throw InvalidParameterException("argument dictionary does not contain the type key");
	}

	euler_type = static_cast<string>(rotation.get_ci("type"));// Warning, will throw


	float e0=0;float e1=0; float e2=0; float e3=0;
	float Omega=0;
	float az  = 0;
	float alt = 0;
	float phi = 0;
	float cxtilt = 0;
	float sxtilt = 0;
	float cytilt = 0;
	float sytilt = 0;
	bool is_quaternion = 0;
	bool is_matrix = 0;

	bool x_mirror;
	float scale;
	// Get these before anything changes so we can apply them again after the rotation is set
	get_scale_and_mirror(scale,x_mirror);
	if (scale == 0) throw UnexpectedBehaviorException("The determinant of the Transform is 0. This is unexpected.");

	string type = Util::str_to_lower(euler_type);
	if (type == "2d") {
		assert_valid_2d();
		az  = 0;
		alt = 0;
		phi = (float)rotation["alpha"] ;
	} else if ( type == "eman" ) {
// 		validate_and_set_type(THREED);
		az  = (float)rotation["az"] ;
		alt = (float)rotation["alt"]  ;
		phi = (float)rotation["phi"] ;
	} else if ( type == "imagic" ) {
// 		validate_and_set_type(THREED);
		az  = (float)rotation["alpha"] ;
		alt = (float)rotation["beta"]  ;
		phi = (float)rotation["gamma"] ;
	} else if ( type == "spider" ) {
// 		validate_and_set_type(THREED);
		az =  (float)rotation["phi"]    + 90.0f;
		alt = (float)rotation["theta"] ;
		phi = (float)rotation["psi"]    - 90.0f;
	} else if ( type == "xyz" ) {
// 		validate_and_set_type(THREED);
		cxtilt = cos( (M_PI/180.0f)*(float)rotation["xtilt"]);
		sxtilt = sin( (M_PI/180.0f)*(float)rotation["xtilt"]);
		cytilt = cos( (M_PI/180.0f)*(float)rotation["ytilt"]);
		sytilt = sin( (M_PI/180.0f)*(float)rotation["ytilt"]);
		az =  (180.0f/M_PI)*atan2(-cytilt*sxtilt,sytilt)   + 90.0f ;
		alt = (180.0f/M_PI)*acos(cytilt*cxtilt)  ;
		phi = (float)rotation["ztilt"] +(180.0f/M_PI)*atan2(sxtilt,cxtilt*sytilt)   - 90.0f ;
	} else if ( type == "mrc" ) {
// 		validate_and_set_type(THREED);
		az  = (float)rotation["phi"]   + 90.0f ;
		alt = (float)rotation["theta"] ;
		phi = (float)rotation["omega"] - 90.0f ;
	} else if ( type == "quaternion" ) {
// 		validate_and_set_type(THREED);
		is_quaternion = 1;
		e0 = (float)rotation["e0"];
		e1 = (float)rotation["e1"];
		e2 = (float)rotation["e2"];
		e3 = (float)rotation["e3"];
	} else if ( type == "spin" ) {
// 		validate_and_set_type(THREED);
		is_quaternion = 1;
		Omega = (float)rotation["Omega"];
		e0 = cos(Omega*M_PI/360.0f);
		e1 = sin(Omega*M_PI/360.0f)* (float)rotation["n1"];
		e2 = sin(Omega*M_PI/360.0f)* (float)rotation["n2"];
		e3 = sin(Omega*M_PI/360.0f)* (float)rotation["n3"];
	} else if ( type == "sgirot" ) {
// 		validate_and_set_type(THREED);
		is_quaternion = 1;
		Omega = (float)rotation["q"] ;
		e0 = cos(Omega*M_PI/360.0f);
		e1 = sin(Omega*M_PI/360.0f)* (float)rotation["n1"];
		e2 = sin(Omega*M_PI/360.0f)* (float)rotation["n2"];
		e3 = sin(Omega*M_PI/360.0f)* (float)rotation["n3"];
	} else if ( type == "matrix" ) {
		is_matrix = 1;
		matrix[0][0] = (float)rotation["m11"];
		matrix[0][1] = (float)rotation["m12"];
		matrix[0][2] = (float)rotation["m13"];
		matrix[1][0] = (float)rotation["m21"];
		matrix[1][1] = (float)rotation["m22"];
		matrix[1][2] = (float)rotation["m23"];
		matrix[2][0] = (float)rotation["m31"];
		matrix[2][1] = (float)rotation["m32"];
		matrix[2][2] = (float)rotation["m33"];
	} else {
// 		transform_type = UNKNOWN;
		throw InvalidStringException(euler_type, "unknown Euler Type");
	}

	float azp  = fmod(az,360.0f)*M_PI/180.0f;
	float altp  = alt*M_PI/180.0f;
	float phip = fmod(phi,360.0f)*M_PI/180.0f;

	if (!is_quaternion && !is_matrix) {
		matrix[0][0] =  cos(phip)*cos(azp) - cos(altp)*sin(azp)*sin(phip);
		matrix[0][1] =  cos(phip)*sin(azp) + cos(altp)*cos(azp)*sin(phip);
		matrix[0][2] =  sin(altp)*sin(phip);
		matrix[1][0] = -sin(phip)*cos(azp) - cos(altp)*sin(azp)*cos(phip);
		matrix[1][1] = -sin(phip)*sin(azp) + cos(altp)*cos(azp)*cos(phip);
		matrix[1][2] =  sin(altp)*cos(phip);
		matrix[2][0] =  sin(altp)*sin(azp);
		matrix[2][1] = -sin(altp)*cos(azp);
		matrix[2][2] =  cos(altp);
	}
	if (is_quaternion){
		matrix[0][0] = e0 * e0 + e1 * e1 - e2 * e2 - e3 * e3;
		matrix[0][1] = 2.0f * (e1 * e2 + e0 * e3);
		matrix[0][2] = 2.0f * (e1 * e3 - e0 * e2);
		matrix[1][0] = 2.0f * (e2 * e1 - e0 * e3);
		matrix[1][1] = e0 * e0 - e1 * e1 + e2 * e2 - e3 * e3;
		matrix[1][2] = 2.0f * (e2 * e3 + e0 * e1);
		matrix[2][0] = 2.0f * (e3 * e1 + e0 * e2);
		matrix[2][1] = 2.0f * (e3 * e2 - e0 * e1);
		matrix[2][2] = e0 * e0 - e1 * e1 - e2 * e2 + e3 * e3;
		// keep in mind matrix[0][2] is M13 gives an e0 e2 piece, etc
	}

	// Apply scale if it existed previously
	if (scale != 1.0f) {
		for(int i=0; i<3; ++i) {
			for(int j=0; j<3; ++j) {
				matrix[i][j] *= scale;
			}
		}
	}

	// Apply post x mirroring if it was applied previouslys
	if ( x_mirror ) {
		for(int j=0; j<3; ++j) {
			matrix[0][j] *= -1.0f;
		}
	}
}

Transform Transform::get_rotation_transform() const
{
	Transform ret(*this);
	ret.set_scale(1.0);
	ret.set_mirror(false);
	ret.set_trans(0,0,0);
	//ret.orthogonalize(); // ?
	return ret;
}

void Transform::set_rotation(const Vec3f & v)
{
	if ( v[0] == 0 && v[1] == 0 && v[2] == 0 )
		throw UnexpectedBehaviorException("Can't set rotation for the null vector");

	Vec3f v1(v);
	v1.normalize();

	float theta = acos(v1[2]); // in radians
	float psi = atan2(v1[1],-v1[0]);
// 	if (v1[1] != 0) psi = asin(v1[1]/sin(theta)); // in radians

	Dict d;
	d["theta"] = (float)EMConsts::rad2deg*theta;
	d["psi"] = (float)EMConsts::rad2deg*psi;
	d["phi"] = (float)0.0;
	d["type"] = "spider";

	set_rotation(d);


}

Transform Transform::negate() const
{
	Transform t(*this);
	for(unsigned int i = 0; i < 3; ++i) {
		for(unsigned int j = 0; j < 4; ++j) {
			t.set(i,j,t[i][j]*-1);
		}
	}
	return t;
}

Transform Transform::get_hflip_transform() const {

	Dict rot = get_rotation("eman");
	rot["alt"] = 180.0f + static_cast<float>(rot["alt"]);
	rot["phi"] = 180.0f - static_cast<float>(rot["phi"]);



	Transform ret(*this); // Is the identity
	ret.set_rotation(rot);

	Vec3f trans = get_trans();
	trans[0] = -trans[0];
	ret.set_trans(trans);

// 	ret.set_mirror(self.get_mirror());

	return ret;
}

Transform Transform::get_vflip_transform() const {

	Dict rot = get_rotation("eman");
	rot["alt"] = 180.0f + static_cast<float>(rot["alt"]);
	rot["phi"] = - static_cast<float>(rot["phi"]);

	Transform ret(*this);
	ret.set_rotation(rot);

	Vec3f trans = get_trans();
	trans[1] = -trans[1];
	ret.set_trans(trans);

	return ret;
}

Dict Transform::get_rotation(const string& euler_type) const
{
	Dict result;

	float max = 1 - ERR_LIMIT;
	float scale;
	bool x_mirror;
	get_scale_and_mirror(scale,x_mirror);
	if (scale == 0) throw UnexpectedBehaviorException("The determinant of the Transform is 0. This is unexpected.");
	float cosalt=matrix[2][2]/scale;
	float x_mirror_scale = (x_mirror ? -1.0f : 1.0f);
	float inv_scale = 1.0f/scale;

	float az=0;
	float alt = 0;
	float phi=0;
	float phiS = 0; // like az   (but in SPIDER ZXZ)
	float psiS =0;  // like phi  (but in SPIDER ZYZ)

	// get alt, az, phi in EMAN convention

	if (cosalt >= 1) {  // that is, alt close to 0 or pi
			alt = 0;
			az=0;
			phi = (float)EMConsts::rad2deg*(float)atan2(x_mirror_scale*matrix[0][1], x_mirror_scale*matrix[0][0]);
	} 
	else if (cosalt <=-1) {
			alt = 180;
			az=0;
			phi=360.0f-(float)EMConsts::rad2deg*(float)atan2(x_mirror_scale*matrix[0][1], x_mirror_scale*matrix[0][0]);
	} else {   // for non exceptional cases:  alt well within 0, pi interval

/*		alt = (float)EMConsts::rad2deg*(float) acos(cosalt);
		az  = 360.0f+(float)EMConsts::rad2deg*(float)atan2(matrix[2][0], -matrix[2][1]);
		phi = 360.0f+(float)EMConsts::rad2deg*(float)atan2(x_mirror_scale*matrix[0][2], matrix[1][2]);// CHANGE ME: This way also works
		phi = 360.0f+(float)EMConsts::rad2deg*(float)atan2(matrix[0][2], matrix[1][2]);  // CHANGE ME: This way works
	//	phi = 360.0f+(float)EMConsts::rad2deg*(float)atan2(x_mirror_scale*matrix[0][2], matrix[1][2]);
*/

		az  = 360.0f+(float)EMConsts::rad2deg*(float)atan2(scale*matrix[2][0], -scale*matrix[2][1]);
//		if (sin(az* (float)EMConsts::deg2rad )* matrix[2][0]< 0) az=az+180.0f;
		alt =     (float)EMConsts::rad2deg* atan(sqrt(matrix[2][0]*matrix[2][0]+matrix[2][1]*matrix[2][1])/fabs(matrix[2][2]) );
		if (matrix[2][2] *scale < 0 ) {
			alt=180.0f-alt;  }
		phi = 360.0f+(float)EMConsts::rad2deg*(float)atan2(x_mirror_scale*matrix[0][2], matrix[1][2]);
	} // ends separate cases: alt close to 0, pi, or neither


// 	az =fmod( az+180.0f,360.0f)-180.0f;
// 	phi=fmod(phi+180.0f,360.0f)-180.0f;

	if (phi > 0) phi = fmod(phi,360.f);
	else phi = 360.f - fmod(fabs(phi),360.f);
	if (phi == 360.f) phi = 0.f;

	if (az > 0 ) az = fmod(az,360.f);
	else az = 360.f - fmod(fabs(az),360.f);
	if (az == 360.f) az = 0.f;

//   get phiS, psiS ; SPIDER
	if (fabs(cosalt) > 1) {  // that is, alt close to 0
		phiS=0;              // used to be if cosalt > max...
		psiS = az+phi;
	}
	else {
		phiS = az   - 90.0f;
		psiS = phi  + 90.0f;
	}
	phiS = fmod((phiS   + 360.0f ), 360.0f) ;
	psiS = fmod((psiS   + 360.0f ), 360.0f) ;

//   do some quaternionic stuff here

	float nphi = (az-phi)/2.0f;
    // The next is also e0
	float cosOover2 = (cos((az+phi)*M_PI/360) * cos(alt*M_PI/360)) ;
	float sinOover2 = sqrt(1 -cosOover2*cosOover2);
	float cosnTheta = sin((az+phi)*M_PI/360) * cos(alt*M_PI/360) / sqrt(1-cosOover2*cosOover2) ;
	float sinnTheta = sqrt(1-cosnTheta*cosnTheta);
	float n1 = sinnTheta*cos(nphi*M_PI/180);
	float n2 = sinnTheta*sin(nphi*M_PI/180);
	float n3 = cosnTheta;
	float xtilt = 0;
	float ytilt = 0;
	float ztilt = 0;


	if (cosOover2<0) {
		cosOover2*=-1; n1 *=-1; n2*=-1; n3*=-1;
	}

	string type = Util::str_to_lower(euler_type);

	result["type"] = type;
	if (type == "2d") {
		assert_valid_2d();
		result["alpha"]  = phi;
	} else if (type == "eman") {
// 		assert_consistent_type(THREED);
		result["az"]  = az;
		result["alt"] = alt;
		result["phi"] = phi;
	} else if (type == "imagic") {
// 		assert_consistent_type(THREED);
		result["alpha"] = az;
		result["beta"] = alt;
		result["gamma"] = phi;
	} else if (type == "spider") {
// 		assert_consistent_type(THREED);
		result["phi"]   = phiS;  // The first Euler like az
		result["theta"] = alt;
		result["psi"]   = psiS;
	} else if (type == "mrc") {
// 		assert_consistent_type(THREED);
		result["phi"]   = phiS;
		result["theta"] = alt;
		result["omega"] = psiS;
	} else if (type == "xyz") {
// 		assert_consistent_type(THREED);
		xtilt = atan2(-sin((M_PI/180.0f)*phiS)*sin((M_PI/180.0f)*alt),cos((M_PI/180.0f)*alt));
		ytilt = asin(  cos((M_PI/180.0f)*phiS)*sin((M_PI/180.0f)*alt));
		ztilt = psiS*M_PI/180.0f - atan2(sin(xtilt), cos(xtilt) *sin(ytilt));

		xtilt=fmod(xtilt*180/M_PI+540.0f,360.0f) -180.0f;
		ztilt=fmod(ztilt*180/M_PI+540.0f,360.0f) -180.0f;

		result["xtilt"]  = xtilt;
		result["ytilt"]  = ytilt*180/M_PI;
		result["ztilt"]  = ztilt;
	} else if (type == "quaternion") {
// 		assert_consistent_type(THREED);
		result["e0"] = cosOover2 ;
		result["e1"] = sinOover2 * n1 ;
		result["e2"] = sinOover2 * n2;
		result["e3"] = sinOover2 * n3;
	} else if (type == "spin") {
// 		assert_consistent_type(THREED);
		result["Omega"] =360.0f* acos(cosOover2)/ M_PI ;
		result["n1"] = n1;
		result["n2"] = n2;
		result["n3"] = n3;
	} else if (type == "sgirot") {
// 		assert_consistent_type(THREED);
		result["q"] = 360.0f*acos(cosOover2)/M_PI ;
		result["n1"] = n1;
		result["n2"] = n2;
		result["n3"] = n3;
	} else if (type == "matrix") {
// 		assert_consistent_type(THREED);
		result["m11"] = x_mirror_scale*matrix[0][0]*inv_scale;
		result["m12"] = x_mirror_scale*matrix[0][1]*inv_scale;
		result["m13"] = x_mirror_scale*matrix[0][2]*inv_scale;
		result["m21"] = matrix[1][0]*inv_scale;
		result["m22"] = matrix[1][1]*inv_scale;
		result["m23"] = matrix[1][2]*inv_scale;
		result["m31"] = matrix[2][0]*inv_scale;
		result["m32"] = matrix[2][1]*inv_scale;
		result["m33"] = matrix[2][2]*inv_scale;
	} else {
		throw InvalidStringException(euler_type, "unknown Euler Type");
	}

	return result;
}

void Transform::set_trans(const float& x, const float& y, const float& z)
{
// 	if ( z != 0.0 ) {
// 		validate_and_set_type(THREED);
// 	}
	bool x_mirror = get_mirror();

	if (x_mirror) matrix[0][3] = -x;
	else matrix[0][3] = x;
	matrix[1][3] = y;
	matrix[2][3] = z;
}

Vec3f Transform::get_trans() const
{
	// No type asserted
	bool x_mirror = get_mirror();
	Vec3f v;
	if (x_mirror) v[0] = -matrix[0][3];
	else v[0] = matrix[0][3];
	v[1] = matrix[1][3];
	v[2] = matrix[2][3];

	Util::apply_precision(v[0],ERR_LIMIT);
	Util::apply_precision(v[1],ERR_LIMIT);
	Util::apply_precision(v[2],ERR_LIMIT);

	return v;
}

Vec2f Transform::get_trans_2d() const
{
	bool x_mirror = get_mirror();
	Vec2f v;
	if (x_mirror) v[0] = -matrix[0][3];
	else v[0] = matrix[0][3];
	v[1] = matrix[1][3];
	return v;
}



Vec3f Transform::get_pre_trans() const
{
	Transform T(*this);
	T.set_trans(0,0,0);
	T.invert();

	Transform soln  = T*(*this);
// 	soln.printme();
	return soln.get_trans();
}

Vec2f Transform::get_pre_trans_2d() const
{
	Transform T(*this);
	T.set_trans(0,0,0);
	T.invert();

	Transform soln  = T*(*this);
// 	soln.printme();
	return soln.get_trans_2d();
}

// Vec2f Transform::get_pre_trans_2d() const
// {
// 	Vec3f v = get_pre_trans();
// 	return Vec2f(v[0],v[1]);
// }


void Transform::set_scale(const float& new_scale) {
	if (new_scale <= 0) {
		throw InvalidValueException(new_scale,"The scale factor in a Transform object must be positive and non zero");
	}
	// Transform = MTSR (Mirroring, Translation, Scaling, Rotate)
	// So changing the scale boils down to this....

	float old_scale = get_scale();

	float n_scale = new_scale;
	Util::apply_precision(n_scale,ERR_LIMIT);

	float corrected_scale = n_scale/old_scale;
	if ( corrected_scale != 1.0 ) {
		for(int i = 0; i < 3;  ++i ) {
			for(int j = 0; j < 3; ++j ) {
				matrix[i][j] *= corrected_scale;
			}
		}
	}
}

float Transform::get_scale() const {
	float determinant = get_determinant();
	if (determinant < 0 ) determinant *= -1;

	float scale = std::pow(determinant,1.0f/3.0f);
	int int_scale = static_cast<int>(scale);
	float scale_residual = scale-static_cast<float>(int_scale);
	if  ( scale_residual < ERR_LIMIT ) { scale = static_cast<float>(int_scale); };

	Util::apply_precision(scale, ERR_LIMIT);

	return scale;
}

void print_matrix(gsl_matrix* M, unsigned int r, unsigned int c, const string& message ) {
	cout << "Message is " << message << endl;
	for ( unsigned int i = 0; i < r; ++i )
	{
		for ( unsigned int j = 0; j < c; ++j )
		{
			cout << gsl_matrix_get(M,i,j) << " ";
		}
		cout << endl;
	}
}

void Transform::orthogonalize()
{
	float scale;
	bool x_mirror;
	get_scale_and_mirror(scale,x_mirror);
	if (scale == 0) throw UnexpectedBehaviorException("The determinant of the Transform is 0. This is unexpected.");
	double inv_scale = 1.0/static_cast<double>(scale);
	double mirror_scale = (x_mirror == true ? -1.0:1.0);

	gsl_matrix * R = gsl_matrix_calloc(3,3);
	for ( unsigned int i = 0; i < 3; ++i )
	{
		for ( unsigned int j = 0; j < 3; ++j )
		{
			if (i == 0 && mirror_scale != 1.0 ) {
				gsl_matrix_set( R, i, j, static_cast<double>(matrix[i][j])*mirror_scale*inv_scale );
			}
			else {
				gsl_matrix_set( R, i, j, static_cast<double>(matrix[i][j])*inv_scale );
			}
		}
	}

	gsl_matrix * V = gsl_matrix_calloc(3,3);
	gsl_vector * S = gsl_vector_calloc(3);
	gsl_vector * work = gsl_vector_calloc(3);
	gsl_linalg_SV_decomp (R, V, S, work); // Now R is U of the SVD R = USV^T

	gsl_matrix * Soln = gsl_matrix_calloc(3,3);
	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, R, V, 0.0, Soln);

	for ( unsigned int i = 0; i < 3; ++i )
	{
		for ( unsigned int j = 0; j < 3; ++j )
		{
			matrix[i][j] = static_cast<float>( gsl_matrix_get(Soln,i,j) );
		}
	}

	// Apply scale if it existed previously
	if (scale != 1.0f) {
		for(int i=0; i<3; ++i) {
			for(int j=0; j<3; ++j) {
				matrix[i][j] *= scale;
			}
		}
	}

	// Apply post x mirroring if it was applied previouslys
	if ( x_mirror ) {
		for(int j=0; j<3; ++j) {
			matrix[0][j] *= -1.0f;
		}
	}

	gsl_matrix_free(V); gsl_matrix_free(R); gsl_matrix_free(Soln);
	gsl_vector_free(S); gsl_vector_free(work);
}

void Transform::set_mirror(const bool x_mirror ) {

	bool old_x_mirror = get_mirror();
	if (old_x_mirror == x_mirror) return; // The user is setting the same value
	else {
		// Toggle the mirroring operation
		for (int j = 0; j < 4; ++j ) {
			matrix[0][j] *= -1;
		}
	}
}

bool Transform::get_mirror() const {
	float determinant = get_determinant();

	bool x_mirror = false;
	if ( determinant < 0 ) x_mirror = true;

	return x_mirror;

}

void Transform::get_scale_and_mirror(float& scale, bool& x_mirror) const {

	float determinant = get_determinant();
	x_mirror = false;
	if ( determinant < 0 ) {
		x_mirror = true;
		determinant *= -1;
	}
 	if (determinant != 1 ) {
		scale = std::pow(determinant,1.0f/3.0f);
		int int_scale = static_cast<int>(scale);
		float scale_residual = scale-static_cast<float>(int_scale);
		if  ( scale_residual < ERR_LIMIT ) { scale = static_cast<float>(int_scale); };
	}
	else scale = 1;

	Util::apply_precision(scale,ERR_LIMIT);
}

float Transform::get_determinant() const
{
	float det = matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[2][1]*matrix[1][2]);
	det -= matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[2][0]*matrix[1][2]);
	det += matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[2][0]*matrix[1][1]);

	Util::apply_precision(det,ERR_LIMIT);

	return det;
}

void Transform::invert() {

	float m00 = matrix[0][0]; float m01=matrix[0][1]; float m02=matrix[0][2];
	float m10 = matrix[1][0]; float m11=matrix[1][1]; float m12=matrix[1][2];
	float m20 = matrix[2][0]; float m21=matrix[2][1]; float m22=matrix[2][2];
	float v0  = matrix[0][3]; float v1 =matrix[1][3]; float v2 =matrix[2][3];

	float cof00 = m11*m22-m12*m21;
	float cof11 = m22*m00-m20*m02;
	float cof22 = m00*m11-m01*m10;
	float cof01 = m10*m22-m20*m12;
	float cof02 = m10*m21-m20*m11;
	float cof12 = m00*m21-m01*m20;
	float cof10 = m01*m22-m02*m21;
	float cof20 = m01*m12-m02*m11;
	float cof21 = m00*m12-m10*m02;

	float det = m00* cof00 + m02* cof02 -m01*cof01;

	matrix[0][0] =   cof00/det;
	matrix[0][1] = - cof10/det;
	matrix[0][2] =   cof20/det;
	matrix[1][0] = - cof01/det;
	matrix[1][1] =   cof11/det;
	matrix[1][2] = - cof21/det;
	matrix[2][0] =   cof02/det;
	matrix[2][1] = - cof12/det;
	matrix[2][2] =   cof22/det;

	matrix[0][3] =  (- cof00*v0 + cof10*v1 - cof20*v2 )/det;
	matrix[1][3] =  (  cof01*v0 - cof11*v1 + cof21*v2 )/det;
	matrix[2][3] =  (- cof02*v0 + cof12*v1 - cof22*v2 )/det;
}

Transform Transform::inverse() const {
	Transform t(*this);
	t.invert();
	return t;
}

void Transform::transpose_inplace() {
	float tempij;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < i; j++) {
			if (i != j) {
				tempij= matrix[i][j];
				matrix[i][j] = matrix[j][i];
				matrix[j][i] = tempij;
			}
		}
	}
}

Transform Transform::transpose() const {
	Transform t(*this);
	t.transpose_inplace();
	return t;
}


Transform EMAN::operator*(const Transform & M2, const Transform & M1)     // YYY
{
	Transform result;
	for (int i=0; i<3; i++) {
		for (int j=0; j<4; j++) {
			result[i][j] = M2[i][0] * M1[0][j] +  M2[i][1] * M1[1][j] + M2[i][2] * M1[2][j];
		}
		result[i][3] += M2[i][3];
	}

	return result;
}

void Transform::assert_valid_2d() const {
	int rotation_error = 0;
	int translation_error = 0;
	if (fabs(matrix[2][0]) > ERR_LIMIT) rotation_error++;
	if (fabs(matrix[2][1]) > ERR_LIMIT) rotation_error++;
	if (fabs(matrix[2][3]) > ERR_LIMIT) translation_error++;
	if (fabs(matrix[0][2]) > ERR_LIMIT) rotation_error++;
	if (fabs(matrix[1][2]) > ERR_LIMIT) rotation_error++;
	if (fabs(matrix[2][2]-1.0) >ERR_LIMIT) rotation_error++; 
	if ( translation_error && rotation_error ) {
		throw UnexpectedBehaviorException("Error, the internal matrix contains 3D rotations and 3D translations. This object can not be considered 2D");
	} else if ( translation_error ) {
		throw UnexpectedBehaviorException("Error, the internal matrix contains a non zero z component for a 3D translation. This object can not be considered 2D");
	}
	else if ( rotation_error ) {
		throw UnexpectedBehaviorException("Error, the internal matrix contains 3D rotations and this object can not be considered 2D");
	}

}



Transform Transform::get_sym(const string & sym_name, int n) const
{
	Symmetry3D* sym = Factory<Symmetry3D>::get(sym_name);
	Transform ret;
	ret = (*this) * sym->get_sym(n);
	delete sym;
	return ret;
}

int Transform::get_nsym(const string & sym_name)
{
	Symmetry3D* sym = Factory<Symmetry3D>::get(sym_name);
	int nsym = sym->get_nsym();
	delete sym;
	return nsym;
}

Transform3D::Transform3D()  //    C1
{
	init();
}

Transform3D::Transform3D( const Transform3D& rhs )
{
    for( int i=0; i < 4; ++i )
    {
        for( int j=0; j < 4; ++j )
	{
	    matrix[i][j] = rhs.matrix[i][j];
	}
    }
}

// C2
Transform3D::Transform3D(const float& az, const float& alt, const float& phi)
{
	init();
	set_rotation(az,alt,phi);
}


//  C3  Usual Constructor: Post Trans, after appying Rot
Transform3D::Transform3D(const float& az, const float& alt, const float& phi, const Vec3f& posttrans )
{
	init(); // This is called in set_rotation
	set_rotation(az,alt,phi);
	set_posttrans(posttrans);
}

Transform3D::Transform3D(const float& m11, const float& m12, const float& m13,
						 const float& m21, const float& m22, const float& m23,
	   const float& m31, const float& m32, const float& m33)
{
	init();
	set_rotation(m11,m12,m13,m21,m22,m23,m31,m32,m33);
}

// C4
Transform3D::Transform3D(EulerType euler_type, const float& a1, const float& a2, const float& a3)
{
 	init();
	set_rotation(euler_type,a1,a2,a3);
}

Transform3D::Transform3D(EulerType euler_type, const float& a1, const float& a2, const float& a3, const float& a4)
{
 	init();
	set_rotation(euler_type,a1,a2,a3,a4);
}


// C5
Transform3D::Transform3D(EulerType euler_type, const Dict& rotation)  //YYY
{
	init();
	set_rotation(euler_type,rotation);
}


// C6   First apply pretrans: Then rotation: Then posttrans

Transform3D::Transform3D(  const Vec3f& pretrans,  const float& az, const float& alt, const float& phi, const Vec3f& posttrans )  //YYY  by default EMAN
{
	init();
	set_pretrans(pretrans);
	set_rotation(az,alt,phi);
	set_posttrans(posttrans);
}




Transform3D::~Transform3D()
{
}



void Transform3D::to_identity()
{
//	for (int i = 0; i < 3; i++) {
//		matrix[i][i] = 1;
//	}

	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			if(i==j) {
				matrix[i][j] = 1;
			}
			else {
				matrix[i][j] = 0;
			}
		}
	}
	post_x_mirror = false;
	set_center(Vec3f(0,0,0));
}



bool Transform3D::is_identity()  // YYY
{
	for (int i=0; i<4; i++) {
		for (int j=0; j<4; j++) {
			if (i==j && matrix[i][j]!=1.0) return 0;
			if (i!=j && matrix[i][j]!=0.0) return 0;
		}
	}
	return 1;
}


void Transform3D::set_center(const Vec3f & center) //YYN
{
	set_pretrans( Vec3f(0,0,0)-center);
	for (int i = 0; i < 3; i++) {
		matrix[i][3]=center[i];
	}
}

//            METHODS
//   Note Transform3Ds are initialized as identities
void Transform3D::init()  // M1
{
	to_identity();
}

//      Set Methods

void Transform3D::set_pretrans(const float& dx, const float& dy, const float& dz) // YYY
{    set_pretrans( Vec3f(dx,dy,dz)); }


void Transform3D::set_pretrans(const float& dx, const float& dy) // YYY
{    set_pretrans( Vec3f(dx,dy,0)); }

void Transform3D::set_pretrans(const Vec2f& pretrans) // YYY
{    set_pretrans( Vec3f(pretrans[0],pretrans[1],0)); }

void Transform3D::set_pretrans(const Vec3f & preT)  // flag=1 means keep the old value of total trans
{
		int flag=0;

//     transFinal = transPost +  Rotation * transPre;
//    This will keep the old value of transPost and change the value of pretrans and the total matrix
    if (flag==0){
		matrix[0][3] = matrix[3][0] + matrix[0][0]*preT[0] + matrix[0][1]*preT[1] + matrix[0][2]*preT[2]  ;
		matrix[1][3] = matrix[3][1] + matrix[1][0]*preT[0] + matrix[1][1]*preT[1] + matrix[1][2]*preT[2]  ;
		matrix[2][3] = matrix[3][2] + matrix[2][0]*preT[0] + matrix[2][1]*preT[1] + matrix[2][2]*preT[2]  ;
	}
//    This will keep the old value of total translation and change the value of posttrans
    if (flag==1){
		matrix[3][0] = matrix[0][3] - (matrix[0][0]*preT[0] + matrix[0][1]*preT[1] + matrix[0][2]*preT[2])  ;
		matrix[3][1] = matrix[1][3] - (matrix[1][0]*preT[0] + matrix[1][1]*preT[1] + matrix[1][2]*preT[2])  ;
		matrix[3][2] = matrix[2][3] - (matrix[2][0]*preT[0] + matrix[2][1]*preT[1] + matrix[2][2]*preT[2])  ;
	}
}


void Transform3D::set_posttrans(const float& dx, const float& dy, const float& dz) // YYY
{    set_posttrans( Vec3f(dx,dy,dz)); }


void Transform3D::set_posttrans(const float& dx, const float& dy) // YYY
{    set_posttrans( Vec3f(dx,dy,0)); }

void Transform3D::set_posttrans(const Vec2f& posttrans) // YYY
{    set_pretrans( Vec3f(posttrans[0],posttrans[1],0)); }

void Transform3D::set_posttrans(const Vec3f & posttrans) // flag=1 means keep the old value of total trans
{
	int flag=0;
    Vec3f preT   = get_pretrans(0) ;
	for (int i = 0; i < 3; i++) {
		matrix[3][i] = posttrans[i];
	}
//     transFinal = transPost +  Rotation * transPre;
//   This will keep the old value of pretrans and change the value of posttrans and the total matrix
	if (flag==0) {
		matrix[0][3] = matrix[3][0] + matrix[0][0]*preT[0] + matrix[0][1]*preT[1] + matrix[0][2]*preT[2]  ;
		matrix[1][3] = matrix[3][1] + matrix[1][0]*preT[0] + matrix[1][1]*preT[1] + matrix[1][2]*preT[2]  ;
		matrix[2][3] = matrix[3][2] + matrix[2][0]*preT[0] + matrix[2][1]*preT[1] + matrix[2][2]*preT[2]  ;
	}
//   This will keep the old value of the total matrix, and c
	if (flag==1) { // Don't do anything
	}
}




void Transform3D::apply_scale(const float& scale)    // YYY
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			matrix[i][j] *= scale;
		}
	}
	for (int j = 0; j < 3; j++) {
		matrix[3][j] *= scale;
	}
}

void Transform3D::orthogonalize()  // YYY
{
	//EulerType EMAN;
	float scale = get_scale() ;
	float inverseScale= 1/scale ;
	apply_scale(inverseScale);
//	Dict angs = get_rotation(EMAN);
//	set_Rotation(EMAN,angs);
}


void Transform3D::transpose()  // YYY
{
	float tempij;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < i; j++) {
			tempij= matrix[i][j];
			matrix[i][j] = matrix[j][i];
			matrix[j][i] = tempij;
		}
	}
}

void Transform3D::set_scale(const float& scale)    // YYY
{
	float OldScale= get_scale();
	float Scale2Apply = scale/OldScale;
	apply_scale(Scale2Apply);
}

float Transform3D::get_mag() const //
{
	EulerType eulertype= SPIN ;
	Dict AA= get_rotation(eulertype);
	return AA["Omega"];
}

Vec3f Transform3D::get_finger() const //
{
	EulerType eulertype= SPIN ;
	Dict AA= get_rotation(eulertype);
	return Vec3f(AA["n1"],AA["n2"],AA["n3"]);
}

Vec3f Transform3D::get_posttrans(int flag) const    //
{
	if (flag==0){
		return Vec3f(matrix[3][0], matrix[3][1], matrix[3][2]);
	}
	// otherwise as if all the translation was post
	return Vec3f(matrix[0][3], matrix[1][3], matrix[2][3]);
}

Vec3f Transform3D::get_total_posttrans() const {
	return get_posttrans(1);
}

Vec3f Transform3D::get_total_pretrans() const {
	return get_pretrans(1);
}


Vec3f Transform3D::get_pretrans(int flag) const    // Fix Me
{
//	The expression is R^T(v_total - v_post);

	Vec3f pretrans;
	Vec3f posttrans(matrix[3][0], matrix[3][1], matrix[3][2]);
	Vec3f tottrans(matrix[0][3], matrix[1][3], matrix[2][3]);
	Vec3f totminuspost;

	totminuspost = tottrans;
	if (flag==0) {
		totminuspost = tottrans-posttrans;
	}

	Transform3D Rinv = inverse();
	for (int i=0; i<3; i++) {
                float ptnow=0;
		for (int j=0; j<3; j++) {
			ptnow +=   Rinv.matrix[i][j]* totminuspost[j] ;
		}
		pretrans.set_value_at(i,ptnow) ;  //
	}
	return pretrans;
}


 Vec3f Transform3D::get_center() const  // YYY
 {
 	return Vec3f();
 }



Vec3f Transform3D::get_matrix3_col(int i) const     // YYY
{
	return Vec3f(matrix[0][i], matrix[1][i], matrix[2][i]);
}


Vec3f Transform3D::get_matrix3_row(int i) const     // YYY
{
	return Vec3f(matrix[i][0], matrix[i][1], matrix[i][2]);
}

Vec3f Transform3D::transform(const Vec3f & v3f) const     // YYY
{
//      This is the transformation of a vector, v by a matrix M
	float x = matrix[0][0] * v3f[0] + matrix[0][1] * v3f[1] + matrix[0][2] * v3f[2] + matrix[0][3] ;
	float y = matrix[1][0] * v3f[0] + matrix[1][1] * v3f[1] + matrix[1][2] * v3f[2] + matrix[1][3] ;
	float z = matrix[2][0] * v3f[0] + matrix[2][1] * v3f[1] + matrix[2][2] * v3f[2] + matrix[2][3] ;
	return Vec3f(x, y, z);
}


Vec3f Transform3D::rotate(const Vec3f & v3f) const     // YYY
{
//      This is the rotation of a vector, v by a matrix M
	float x = matrix[0][0] * v3f[0] + matrix[0][1] * v3f[1] + matrix[0][2] * v3f[2]  ;
	float y = matrix[1][0] * v3f[0] + matrix[1][1] * v3f[1] + matrix[1][2] * v3f[2]  ;
	float z = matrix[2][0] * v3f[0] + matrix[2][1] * v3f[1] + matrix[2][2] * v3f[2]  ;
	return Vec3f(x, y, z);
}


Transform3D EMAN::operator*(const Transform3D & M2, const Transform3D & M1)     // YYY
{
//       This is the  left multiplication of a matrix M1 by a matrix M2; that is M2*M1
//       It returns a new matrix
	Transform3D resultant;
	for (int i=0; i<3; i++) {
		for (int j=0; j<4; j++) {
			resultant[i][j] = M2[i][0] * M1[0][j] +  M2[i][1] * M1[1][j] + M2[i][2] * M1[2][j];
		}
		resultant[i][3] += M2[i][3];  // add on the new translation (not included above)
	}

	for (int j=0; j<3; j++) {
		resultant[3][j] = M2[3][j];
	}

	return resultant; // This will have the post_trans of M2
}

/*             Here starts the pure rotation stuff */

/**  A rotation is given by

EMAN
  | cos phi   sin phi    0 |  |  1       0     0      | |  cos az  sin az   0 |
  |-sin phi   cos phi    0 |  |  0   cos alt  sin alt | | -sin az  cos az   0 |
  |   0          0       1 |  |  0  -sin alt  cos alt | |     0       0     1 |

---------------------------------------------------------------------------

SPIDER, FREEALIGN  (th == theta)
| cos psi   sin psi    0 |  |  cos th  0   -sin th   | |  cos phi  sin phi   0 |
|-sin psi   cos psi    0 |  |  0       1       0     | | -sin phi  cos phi   0 |
|   0          0       1 |  |  sin th  0    cos th   | |     0        0      1 |


Now this middle matrix is equal to

                 | 0 -1 0|  |1     0    0     | | 0  1  0 |
                 | 1  0 0|  |0  cos th sin th | |-1  0  0 |
                 | 0  0 1|  |0 -sin th cos th | | 0  0  1 |


 So we have

  | sin psi  -cos psi    0 |  |  1       0     0    | | -sin phi  cos phi   0 |
  | cos psi   sin psi    0 |  |  0   cos th  sin th | | -cos phi -sin phi   0 |
  |   0          0       1 |  |  0  -sin th  cos th | |     0       0     1 |


        so az = phi_SPIDER + pi/2
          phi = psi        - pi/2

---------------------------------------------------------------------------

MRC  th=theta; om=omega ;

dwoolford says - this is wrong, the derivation of phi is the negative of the true result

| cos om   sin om    0 |  |  cos th  0   -sin th   | |  cos phi  sin phi   0 |
|-sin om   cos om    0 |  |  0       1       0     | | -sin phi  cos phi   0 |
|   0        0       1 |  |  sin th  0    cos th   | |     0        0      1 |

        so az = phi     + pi/2
          alt = theta
          phi = omega   - pi/2

---------------------------------------------------------------------------
For the quaternion type operations, we can start with

R = (1-nhat nhat) cos(Omega) - sin(Omega)nhat cross + nhat nhat
Notice that this is a clockwise rotation( the xy component, for nhat=zhat,
 is calculated as - sin(Omega) xhat dot zhat cross yhat= sin(Omega): this is the
 correct sign for clockwise rotations).
Now we develop

R =  cos(Omega) one + nhat nhat (1-cos(Omega)) - sin(Omega) nhat cross
  = (cos^2(Omega/2) - sin^2(Omega/2)) one  + 2 ((sin(Omega/2)nhat ) ((sin(Omega/2)nhat )
                                    - 2 cos(Omega/2) ((sin(Omega/2)nhat )  cross
  = (e0^2 - evec^2) one  + 2 (evec evec )  - 2 e0 evec  cross

  e0 = cos(Omega/2)
  vec{e} = sin(Omega/2) nhat


SGIrot is the same as SPIN (see paper)
The update of rotations for quaternions is very easy.


*/

void Transform3D::set_rotation(const float& az, const float& alt, const float& phi )
{
	EulerType euler_type=EMAN;
	Dict rot;
	rot["az"]  = az;
	rot["alt"] = alt;
	rot["phi"] = phi;
	set_rotation(euler_type, rot);
}

// This is where it all happens;
void Transform3D::set_rotation(EulerType euler_type, const float& a1, const float& a2, const float& a3)
{
	init();
	Dict rot;
	switch(euler_type) {
		case EMAN:
			rot["az"]  = a1;
			rot["alt"] = a2;
			rot["phi"] = a3;
			break;
		case SPIDER:
			rot["phi"]   = a1;
			rot["theta"] = a2;
			rot["psi"]   = a3;
			break;
		case IMAGIC:
			rot["alpha"]   = a1;
			rot["beta"] = a2;
			rot["gamma"]   = a3;
			break;
		case MRC:
			rot["phi"]   = a1;
			rot["theta"] = a2;
			rot["omega"]   = a3;
			break;
		case XYZ:
			rot["xtilt"]   = a1;
			rot["ytilt"] = a2;
			rot["ztilt"]   = a3;
			break;
		default:
		throw InvalidValueException(euler_type, "cannot instantiate this Euler Type");
  	}  // ends switch euler_type
	set_rotation(euler_type, rot);
}

// This is where it all happens;
void Transform3D::set_rotation(EulerType euler_type, const float& a1, const float& a2, const float& a3, const float& a4)
{
	init();
	Dict rot;
	switch(euler_type) {
		case QUATERNION:
			rot["e0"]  = a1;
			rot["e1"] = a2;
			rot["e2"] = a3;
			rot["e3"] = a4;
			break;
		case SGIROT:
			rot["q"]  = a1;
			rot["n1"] = a2;
			rot["n2"] = a3;
			rot["n3"] = a4;
		case SPIN:
			rot["Omega"]  = a1;
			rot["n1"] = a2;
			rot["n2"] = a3;
			rot["n3"] = a4;
			break;
		default:
			throw InvalidValueException(euler_type, "cannot instantiate this Euler Type");
	}  // ends switch euler_type
	set_rotation(euler_type, rot);
}

void Transform3D::set_rotation(EulerType euler_type, const Dict& rotation)
{
	float e0  = 0;float e1=0; float e2=0; float e3=0;
	float Omega=0;
	float az  = 0;
	float alt = 0;
	float phi = 0;
	float cxtilt = 0;
	float sxtilt = 0;
	float cytilt = 0;
	float sytilt = 0;
	bool is_quaternion = 0;
	bool is_matrix = 0;

	switch(euler_type) {
	case EMAN:
		az  = (float)rotation["az"] ;
		alt = (float)rotation["alt"]  ;
		phi = (float)rotation["phi"] ;
		break;
	case IMAGIC:
		az  = (float)rotation["alpha"] ;
		alt = (float)rotation["beta"]  ;
		phi = (float)rotation["gamma"] ;
		break;

	case SPIDER:
		az =  (float)rotation["phi"]    + 90.0f;
		alt = (float)rotation["theta"] ;
		phi = (float)rotation["psi"]    - 90.0f;
		break;

	case XYZ:
		cxtilt = cos( (M_PI/180.0f)*(float)rotation["xtilt"]);
		sxtilt = sin( (M_PI/180.0f)*(float)rotation["xtilt"]);
		cytilt = cos( (M_PI/180.0f)*(float)rotation["ytilt"]);
		sytilt = sin( (M_PI/180.0f)*(float)rotation["ytilt"]);
		az =  (180.0f/M_PI)*atan2(-cytilt*sxtilt,sytilt)   + 90.0f ;
		alt = (180.0f/M_PI)*acos(cytilt*cxtilt)  ;
		phi = (float)rotation["ztilt"] +(180.0f/M_PI)*atan2(sxtilt,cxtilt*sytilt)   - 90.0f ;
		break;

	case MRC:
		az  = (float)rotation["phi"]   + 90.0f ;
		alt = (float)rotation["theta"] ;
		phi = (float)rotation["omega"] - 90.0f ;
		break;

	case QUATERNION:
		is_quaternion = 1;
		e0 = (float)rotation["e0"];
		e1 = (float)rotation["e1"];
		e2 = (float)rotation["e2"];
		e3 = (float)rotation["e3"];
		break;

	case SPIN:
		is_quaternion = 1;
		Omega = (float)rotation["Omega"];
		e0 = cos(Omega*M_PI/360.0f);
		e1 = sin(Omega*M_PI/360.0f)* (float)rotation["n1"];
		e2 = sin(Omega*M_PI/360.0f)* (float)rotation["n2"];
		e3 = sin(Omega*M_PI/360.0f)* (float)rotation["n3"];
		break;

	case SGIROT:
		is_quaternion = 1;
		Omega = (float)rotation["q"]  ;
		e0 = cos(Omega*M_PI/360.0f);
		e1 = sin(Omega*M_PI/360.0f)* (float)rotation["n1"];
		e2 = sin(Omega*M_PI/360.0f)* (float)rotation["n2"];
		e3 = sin(Omega*M_PI/360.0f)* (float)rotation["n3"];
		break;

	case MATRIX:
		is_matrix = 1;
		matrix[0][0] = (float)rotation["m11"]  ;
		matrix[0][1] = (float)rotation["m12"]  ;
		matrix[0][2] = (float)rotation["m13"]  ;
		matrix[1][0] = (float)rotation["m21"]  ;
		matrix[1][1] = (float)rotation["m22"]  ;
		matrix[1][2] = (float)rotation["m23"]  ;
		matrix[2][0] = (float)rotation["m31"]  ;
		matrix[2][1] = (float)rotation["m32"]  ;
		matrix[2][2] = (float)rotation["m33"]  ;
		break;

	default:
		throw InvalidValueException(euler_type, "unknown Euler Type");
	}  // ends switch euler_type


	Vec3f postT  = get_posttrans( ) ;
	Vec3f preT   = get_pretrans( ) ;


	float azp  = fmod(az,360.0f)*M_PI/180.0f;
	float altp  = alt*M_PI/180.0f;
	float phip = fmod(phi,360.0f)*M_PI/180.0f;

	if (!is_quaternion && !is_matrix) {
		matrix[0][0] =  cos(phip)*cos(azp) - cos(altp)*sin(azp)*sin(phip);
		matrix[0][1] =  cos(phip)*sin(azp) + cos(altp)*cos(azp)*sin(phip);
		matrix[0][2] =  sin(altp)*sin(phip);
		matrix[1][0] = -sin(phip)*cos(azp) - cos(altp)*sin(azp)*cos(phip);
		matrix[1][1] = -sin(phip)*sin(azp) + cos(altp)*cos(azp)*cos(phip);
		matrix[1][2] =  sin(altp)*cos(phip);
		matrix[2][0] =  sin(altp)*sin(azp);
		matrix[2][1] = -sin(altp)*cos(azp);
		matrix[2][2] =  cos(altp);
	}
	if (is_quaternion){
		matrix[0][0] = e0 * e0 + e1 * e1 - e2 * e2 - e3 * e3;
		matrix[0][1] = 2.0f * (e1 * e2 + e0 * e3);
		matrix[0][2] = 2.0f * (e1 * e3 - e0 * e2);
		matrix[1][0] = 2.0f * (e2 * e1 - e0 * e3);
		matrix[1][1] = e0 * e0 - e1 * e1 + e2 * e2 - e3 * e3;
		matrix[1][2] = 2.0f * (e2 * e3 + e0 * e1);
		matrix[2][0] = 2.0f * (e3 * e1 + e0 * e2);
		matrix[2][1] = 2.0f * (e3 * e2 - e0 * e1);
		matrix[2][2] = e0 * e0 - e1 * e1 - e2 * e2 + e3 * e3;
		// keep in mind matrix[0][2] is M13 gives an e0 e2 piece, etc
	}
	// Now do post and pretrans: vfinal = vpost + R vpre;

	matrix[0][3] = postT[0] + matrix[0][0]*preT[0] + matrix[0][1]*preT[1] + matrix[0][2]*preT[2]  ;
	matrix[1][3] = postT[1] + matrix[1][0]*preT[0] + matrix[1][1]*preT[1] + matrix[1][2]*preT[2]  ;
	matrix[2][3] = postT[2] + matrix[2][0]*preT[0] + matrix[2][1]*preT[1] + matrix[2][2]*preT[2]  ;
}


void Transform3D::set_rotation(const float& m11, const float& m12, const float& m13,
							   const float& m21, const float& m22, const float& m23,
		  const float& m31, const float& m32, const float& m33)
{
	EulerType euler_type = MATRIX;
	Dict rot;
	rot["m11"]  = m11;
	rot["m12"]  = m12;
	rot["m13"]  = m13;
	rot["m21"]  = m21;
	rot["m22"]  = m22;
	rot["m23"]  = m23;
	rot["m31"]  = m31;
	rot["m32"]  = m32;
	rot["m33"]  = m33;
	set_rotation(euler_type, rot);  // Or should it be &rot ?
}

void Transform3D::set_rotation(const Vec3f & eahat, const Vec3f & ebhat,
                                    const Vec3f & eAhat, const Vec3f & eBhat)
{// this rotation rotates unit vectors a,b into A,B;
//    The program assumes a dot b must equal A dot B
	Vec3f eahatcp(eahat);
	Vec3f ebhatcp(ebhat);
	Vec3f eAhatcp(eAhat);
	Vec3f eBhatcp(eBhat);

	eahatcp.normalize();
	ebhatcp.normalize();
	eAhatcp.normalize();
	eBhatcp.normalize();

	Vec3f aMinusA(eahatcp);
	aMinusA  -= eAhatcp;
	Vec3f bMinusB(ebhatcp);
	bMinusB  -= eBhatcp;


	Vec3f  nhat;
	float aAlength = aMinusA.length();
	float bBlength = bMinusB.length();
	if (aAlength==0){
		nhat=eahatcp;
	}else if (bBlength==0){
		nhat=ebhatcp;
	}else{
		nhat= aMinusA.cross(bMinusB);
		nhat.normalize();
	}

//		printf("nhat=%f,%f,%f \n",nhat[0],nhat[1],nhat[2]);

	Vec3f neahat  = eahatcp.cross(nhat);
	Vec3f nebhat  = ebhatcp.cross(nhat);
	Vec3f neAhat  = eAhatcp.cross(nhat);
	Vec3f neBhat  = eBhatcp.cross(nhat);

	float cosOmegaA = (neahat.dot(neAhat))  / (neahat.dot(neahat));
//	float cosOmegaB = (nebhat.dot(neBhat))  / (nebhat.dot(nebhat));
	float sinOmegaA = (neahat.dot(eAhatcp)) / (neahat.dot(neahat));
//	printf("cosOmegaA=%f \n",cosOmegaA); 	printf("sinOmegaA=%f \n",sinOmegaA);

	float OmegaA = atan2(sinOmegaA,cosOmegaA);
//	printf("OmegaA=%f \n",OmegaA*180/M_PI);

	EulerType euler_type=SPIN;
	Dict rotation;
	rotation["n1"]= nhat[0];
	rotation["n2"]= nhat[1];
	rotation["n3"]= nhat[2];
	rotation["Omega"] =OmegaA*180.0/M_PI;
	set_rotation(euler_type,  rotation);
}


float Transform3D::get_scale() const     // YYY
{
	// Assumes uniform scaling, calculation uses Z only.
	float scale =0;
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			scale = scale + matrix[i][j]*matrix[i][j];
		}
	}

	return sqrt(scale/3);
}



Dict Transform3D::get_rotation(EulerType euler_type) const
{
	Dict result;

	float max = 1 - ERR_LIMIT;
	float sca=get_scale();
	float cosalt=matrix[2][2]/sca;


	float az=0;
	float alt = 0;
	float phi=0;
	float phiS = 0; // like az   (but in SPIDER ZXZ)
	float psiS =0;  // like phi  (but in SPIDER ZYZ)


// get alt, az, phi;  EMAN

	if (cosalt > max) {  // that is, alt close to 0
		alt = 0;
		az=0;
		phi = (float)EMConsts::rad2deg*(float)atan2(matrix[0][1], matrix[0][0]);
	}
	else if (cosalt < -max) { // alt close to pi
		alt = 180;
		az=0;
		phi=360.0f-(float)EMConsts::rad2deg*(float)atan2(matrix[0][1], matrix[0][0]);
	}
	else {
		alt = (float)EMConsts::rad2deg*(float) acos(cosalt);
		az  = 360.0f+(float)EMConsts::rad2deg*(float)atan2(matrix[2][0], -matrix[2][1]);
		phi = 360.0f+(float)EMConsts::rad2deg*(float)atan2(matrix[0][2], matrix[1][2]);
	}
	az=fmod(az+180.0f,360.0f)-180.0f;
	phi=fmod(phi+180.0f,360.0f)-180.0f;

//   get phiS, psiS ; SPIDER
	if (fabs(cosalt) > max) {  // that is, alt close to 0
		phiS=0;
		psiS = az+phi;
	}
	else {
		phiS = az   - 90.0f;
		psiS = phi  + 90.0f;
	}
	phiS = fmod((phiS   + 360.0f ), 360.0f) ;
	psiS = fmod((psiS   + 360.0f ), 360.0f) ;

//   do some quaternionic stuff here

	float nphi = (az-phi)/2.0f;
    // The next is also e0
	float cosOover2 = (cos((az+phi)*M_PI/360) * cos(alt*M_PI/360)) ;
	float sinOover2 = sqrt(1 -cosOover2*cosOover2);
	float cosnTheta = sin((az+phi)*M_PI/360) * cos(alt*M_PI/360) / sqrt(1-cosOover2*cosOover2) ;
	float sinnTheta = sqrt(1-cosnTheta*cosnTheta);
	float n1 = sinnTheta*cos(nphi*M_PI/180);
	float n2 = sinnTheta*sin(nphi*M_PI/180);
	float n3 = cosnTheta;
        float xtilt = 0;
        float ytilt = 0;
        float ztilt = 0;


	if (cosOover2<0) {
		cosOover2*=-1; n1 *=-1; n2*=-1; n3*=-1;
	}


	switch (euler_type) {
	case EMAN:
		result["az"]  = az;
		result["alt"] = alt;
		result["phi"] = phi;
		break;

	case IMAGIC:
		result["alpha"] = az;
		result["beta"] = alt;
		result["gamma"] = phi;
		break;

	case SPIDER:
		result["phi"]   = phiS;  // The first Euler like az
		result["theta"] = alt;
		result["psi"]   = psiS;
		break;

	case MRC:
		result["phi"]   = phiS;
		result["theta"] = alt;
		result["omega"] = psiS;
		break;

	case XYZ:
	        xtilt = atan2(-sin((M_PI/180.0f)*phiS)*sin((M_PI/180.0f)*alt),cos((M_PI/180.0f)*alt));
	        ytilt = asin(  cos((M_PI/180.0f)*phiS)*sin((M_PI/180.0f)*alt));
	        ztilt = psiS*M_PI/180.0f - atan2(sin(xtilt), cos(xtilt) *sin(ytilt));

		xtilt=fmod(xtilt*180/M_PI+540.0f,360.0f) -180.0f;
		ztilt=fmod(ztilt*180/M_PI+540.0f,360.0f) -180.0f;

		result["xtilt"]  = xtilt;
		result["ytilt"]  = ytilt*180/M_PI;
		result["ztilt"]  = ztilt;
		break;

	case QUATERNION:
		result["e0"] = cosOover2 ;
		result["e1"] = sinOover2 * n1 ;
		result["e2"] = sinOover2 * n2;
		result["e3"] = sinOover2 * n3;
		break;

	case SPIN:
		result["Omega"] =360.0f* acos(cosOover2)/ M_PI ;
		result["n1"] = n1;
		result["n2"] = n2;
		result["n3"] = n3;
		break;

	case SGIROT:
		result["q"] = 360.0f*acos(cosOover2)/M_PI ;
		result["n1"] = n1;
		result["n2"] = n2;
		result["n3"] = n3;
		break;

	case MATRIX:
		result["m11"] = matrix[0][0] ;
		result["m12"] = matrix[0][1] ;
		result["m13"] = matrix[0][2] ;
		result["m21"] = matrix[1][0] ;
		result["m22"] = matrix[1][1] ;
		result["m23"] = matrix[1][2] ;
		result["m31"] = matrix[2][0] ;
		result["m32"] = matrix[2][1] ;
		result["m33"] = matrix[2][2] ;
		break;

	default:
		throw InvalidValueException(euler_type, "unknown Euler Type");
	}

	return result;
}

Transform3D Transform3D::inverseUsingAngs() const    //   YYN need to test it for sure
{
	// First Find the scale
	EulerType eE=EMAN;


	float scale   = get_scale();
	Vec3f preT   = get_pretrans( ) ;
	Vec3f postT   = get_posttrans( ) ;
	Dict angs     = get_rotation(eE);
	Dict invAngs  ;

	invAngs["phi"]   = 180.0f - (float) angs["az"] ;
	invAngs["az"]    = 180.0f - (float) angs["phi"] ;
	invAngs["alt"]   = angs["alt"] ;

//    The inverse result
//
//           Z_phi   X_alt     Z_az
//                 is
//       Z_{pi-az}   X_alt  Z_{pi-phi}
//      The reason for the extra pi's, is because one would like to keep alt positive

	float inverseScale= 1/scale ;

	Transform3D invM;

	invM.set_rotation(EMAN, invAngs);
	invM.apply_scale(inverseScale);
	invM.set_pretrans(-postT );
	invM.set_posttrans(-preT );


	return invM;

}

Transform3D Transform3D::inverse() const    //   YYN need to test it for sure
{
	// This assumes the matrix is 4 by 4 and the last row reads [0 0 0 1]

	float m00 = matrix[0][0]; float m01=matrix[0][1]; float m02=matrix[0][2];
	float m10 = matrix[1][0]; float m11=matrix[1][1]; float m12=matrix[1][2];
	float m20 = matrix[2][0]; float m21=matrix[2][1]; float m22=matrix[2][2];
	float v0  = matrix[0][3]; float v1 =matrix[1][3]; float v2 =matrix[2][3];

    float cof00 = m11*m22-m12*m21;
    float cof11 = m22*m00-m20*m02;
    float cof22 = m00*m11-m01*m10;
    float cof01 = m10*m22-m20*m12;
    float cof02 = m10*m21-m20*m11;
    float cof12 = m00*m21-m01*m20;
    float cof10 = m01*m22-m02*m21;
    float cof20 = m01*m12-m02*m11;
    float cof21 = m00*m12-m10*m02;

    float Det = m00* cof00 + m02* cof02 -m01*cof01;

    Transform3D invM;

    invM.matrix[0][0] =   cof00/Det;
    invM.matrix[0][1] = - cof10/Det;
    invM.matrix[0][2] =   cof20/Det;
    invM.matrix[1][0] = - cof01/Det;
    invM.matrix[1][1] =   cof11/Det;
    invM.matrix[1][2] = - cof21/Det;
    invM.matrix[2][0] =   cof02/Det;
    invM.matrix[2][1] = - cof12/Det;
    invM.matrix[2][2] =   cof22/Det;

    invM.matrix[0][3] =  (- cof00*v0 + cof10*v1 - cof20*v2 )/Det;
    invM.matrix[1][3] =  (  cof01*v0 - cof11*v1 + cof21*v2 )/Det;
    invM.matrix[2][3] =  (- cof02*v0 + cof12*v1 - cof22*v2 )/Det;

	Vec3f postT   = get_posttrans( ) ;
	Vec3f invMpre   = - postT;
	Vec3f invMpost   ;
	for ( int i = 0; i < 3; i++) {
		invMpost[i] = invM.matrix[i][3];
		for ( int j = 0; j < 3; j++) {
			invMpost[i] += - invM.matrix[i][j]*invMpre[j];
		}
		invM.matrix[3][i] = invMpost[i];
	}

	return invM;
}



// Symmetry Stuff

Transform3D Transform3D::get_sym(const string & symname, int n) const
{
	int nsym = get_nsym(symname);

// 	Transform3D invalid;
// 	invalid.set_rotation( -0.1f, -0.1f, -0.1f);

	// see www.math.utah.edu/~alfeld/math/polyhedra/polyhedra.html for pictures
	// By default we will put largest symmetry along z-axis.

	// Each Platonic Solid has 2E symmetry elements.


	// An icosahedron has   m=5, n=3, F=20 E=30=nF/2, V=12=nF/m,since vertices shared by 5 triangles;
	// It is composed of 20 triangles. E=3*20/2;


	// An dodecahedron has m=3, n=5   F=12 E=30  V=20
	// It is composed of 12 pentagons. E=5*12/2;   V= 5*12/3, since vertices shared by 3 pentagons;



    // The ICOS symmetry group has the face along z-axis

	float lvl0=0;                             //  there is one pentagon on top; five-fold along z
	float lvl1= 63.4349f; // that is atan(2)  // there are 5 pentagons with centers at this height (angle)
	float lvl2=116.5651f; //that is 180-lvl1  // there are 5 pentagons with centers at this height (angle)
	float lvl3=180.f;                           // there is one pentagon on the bottom
             // Notice that 63.439 is the angle between two faces of the dual object

	static double ICOS[180] = { // This is with a pentagon normal to z
		  0,lvl0,0,    0,lvl0,288,   0,lvl0,216,   0,lvl0,144,  0,lvl0,72,
		  0,lvl1,36,   0,lvl1,324,   0,lvl1,252,   0,lvl1,180,  0,lvl1,108,
		 72,lvl1,36,  72,lvl1,324,  72,lvl1,252,  72,lvl1,180,  72,lvl1,108,
		144,lvl1,36, 144,lvl1,324, 144,lvl1,252, 144,lvl1,180, 144,lvl1,108,
		216,lvl1,36, 216,lvl1,324, 216,lvl1,252, 216,lvl1,180, 216,lvl1,108,
		288,lvl1,36, 288,lvl1,324, 288,lvl1,252, 288,lvl1,180, 288,lvl1,108,
		 36,lvl2,0,   36,lvl2,288,  36,lvl2,216,  36,lvl2,144,  36,lvl2,72,
		108,lvl2,0,  108,lvl2,288, 108,lvl2,216, 108,lvl2,144, 108,lvl2,72,
		180,lvl2,0,  180,lvl2,288, 180,lvl2,216, 180,lvl2,144, 180,lvl2,72,
		252,lvl2,0,  252,lvl2,288, 252,lvl2,216, 252,lvl2,144, 252,lvl2,72,
		324,lvl2,0,  324,lvl2,288, 324,lvl2,216, 324,lvl2,144, 324,lvl2,72,
   		  0,lvl3,0,    0,lvl3,288,   0,lvl3,216,   0,lvl3,144,   0,lvl3,72
	};


	// A cube has   m=3, n=4, F=6 E=12=nF/2, V=8=nF/m,since vertices shared by 3 squares;
	// It is composed of 6 squares.


	// An octahedron has   m=4, n=3, F=8 E=12=nF/2, V=6=nF/m,since vertices shared by 4 triangles;
	// It is composed of 8 triangles.

    // We have placed the OCT symmetry group with a face along the z-axis
        lvl0=0;
	lvl1=90;
	lvl2=180;

	static float OCT[72] = {// This is with a face of a cube along z
		      0,lvl0,0,   0,lvl0,90,    0,lvl0,180,    0,lvl0,270,
		      0,lvl1,0,   0,lvl1,90,    0,lvl1,180,    0,lvl1,270,
		     90,lvl1,0,  90,lvl1,90,   90,lvl1,180,   90,lvl1,270,
		    180,lvl1,0, 180,lvl1,90,  180,lvl1,180,  180,lvl1,270,
		    270,lvl1,0, 270,lvl1,90,  270,lvl1,180,  270,lvl1,270,
		      0,lvl2,0,   0,lvl2,90,    0,lvl2,180,    0,lvl2,270
	};
	// B^4=A^3=1;  BABA=1; implies   AA=BAB, ABA=B^3 , AB^2A = BBBABBB and
	//   20 words with at most a single A
    //   1 B BB BBB A  BA AB BBA BAB ABB BBBA BBAB BABB ABBB BBBAB BBABB BABBB
    //                        BBBABB BBABBB BBBABBB
     // also     ABBBA is distinct yields 4 more words
     //    ABBBA   BABBBA BBABBBA BBBABBBA
     // for a total of 24 words
     // Note A BBB A BBB A  reduces to BBABB
     //  and  B A BBB A is the same as A BBB A BBB etc.

    // The TET symmetry group has a face along the z-axis
    // It has n=m=3; F=4, E=6=nF/2, V=4=nF/m
        lvl0=0;         // There is a face along z
	lvl1=109.4712f;  //  that is acos(-1/3)  // There  are 3 faces at this angle

	static float TET[36] = {// This is with the face along z
	      0,lvl0,0,   0,lvl0,120,    0,lvl0,240,
	      0,lvl1,60,   0,lvl1,180,    0,lvl1,300,
	    120,lvl1,60, 120,lvl1,180,  120,lvl1,300,
	    240,lvl1,60, 240,lvl1,180,  240,lvl1,300
	};
	// B^3=A^3=1;  BABA=1; implies   A^2=BAB, ABA=B^2 , AB^2A = B^2AB^2 and
	//   12 words with at most a single A
    //   1 B BB  A  BA AB BBA BAB ABB BBAB BABB BBABB
    // at most one A is necessary

	Transform3D ret;
	SymType type = get_sym_type(symname);

	switch (type) {
	case CSYM:
		ret.set_rotation( n * 360.0f / nsym, 0, 0);
		break;
	case DSYM:
		if (n >= nsym / 2) {
			ret.set_rotation((n - nsym/2) * 360.0f / (nsym / 2),180.0f, 0);
		}
		else {
			ret.set_rotation( n * 360.0f / (nsym / 2),0, 0);
		}
		break;
	case ICOS_SYM:
		ret.set_rotation((float)ICOS[n * 3 ],
				 (float)ICOS[n * 3 + 1],
				 (float)ICOS[n * 3 + 2] );
		break;
	case OCT_SYM:
		ret.set_rotation((float)OCT[n * 3],
				 (float)OCT[n * 3 + 1],
				 (float)OCT[n * 3 + 2] );
		break;
	case TET_SYM:
		ret.set_rotation((float)TET[n * 3 ],
				 (float)TET[n * 3 + 1] ,
				 (float)TET[n * 3 + 2] );
		break;
	case ISYM:
		ret.set_rotation(0, 0, 0);
		break;
	default:
		throw InvalidValueException(type, symname);
	}

	ret = (*this) * ret;

	return ret;
}

int Transform3D::get_nsym(const string & name)
{
	string symname = name;

	for (size_t i = 0; i < name.size(); i++) {
		if (isalpha(name[i])) {
			symname[i] = (char)tolower(name[i]);
		}
	}

	SymType type = get_sym_type(symname);
	int nsym = 0;

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
	case TET_SYM:
		nsym = 12;
		break;
	case ISYM:
		nsym = 1;
		break;
	case UNKNOWN_SYM:
	default:
		throw InvalidValueException(type, name);
	}
	return nsym;
}



Transform3D::SymType Transform3D::get_sym_type(const string & name)
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
	else if (name == "tet") {
		t = TET_SYM;
	}
	else if (name == "i" || name == "") {
		t = ISYM;
	}
	return t;
}

vector<Transform3D*>
Transform3D::angles2tfvec(EulerType eulertype, const vector<float> ang) {
	int nangles = ang.size() / 3;
	vector<Transform3D*> tfvec;
	for (int i = 0; i < nangles; i++) {
		tfvec.push_back(new Transform3D(eulertype,ang[3*i],ang[3*i+1],ang[3*i+2]));
	}
	return tfvec;
}




/* vim: set ts=4 noet: */
