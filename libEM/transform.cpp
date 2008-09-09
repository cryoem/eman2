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
using namespace EMAN;

#ifdef WIN32
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif

#include <algorithm> // for std::transform

				 
const float Transform3D::ERR_LIMIT = 0.000001f;

const string CSym::NAME = "c";
const string DSym::NAME = "d";
const string HSym::NAME = "h";
const string TetrahedralSym::NAME = "tet";
const string OctahedralSym::NAME = "oct";
const string IcosahedralSym::NAME = "icos";
const string EmanOrientationGenerator::NAME = "eman";
const string SaffOrientationGenerator::NAME = "saff";
const string EvenOrientationGenerator::NAME = "even";
const string RandomOrientationGenerator::NAME = "rand";
const string OptimumOrientationGenerator::NAME = "opt";

const float Transform::ERR_LIMIT = 0.000001f;

 Transform::Transform() : transform_type(UNKNOWN)
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
		transform_type = that.transform_type;
	}	
	return *this;
}

Transform::Transform(const Dict& d) : transform_type(UNKNOWN) {
	to_identity();
	set_params(d);
}


void Transform::to_identity()
{
	transform_type = UNKNOWN;
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

void Transform::set_params(const Dict& d) {
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
		if ( (e.get_type() != EMObjectTypes::BOOL ) && (e.get_type() != EMObjectTypes::INT ) && (e.get_type() != EMObjectTypes::UNSIGNEDINT ) )
			throw InvalidParameterException("Error, mirror must be a bool or an int");
		
		bool mirror = static_cast<bool>(e);
		set_mirror(mirror);
	}

}


Dict Transform::get_params(const string& euler_type) {
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

Dict Transform::get_params_inverse(const string& euler_type) {
	Dict params = get_rotation(euler_type);
	
	Vec3f v = get_pre_trans();
	params["tx"] = v[0]; params["ty"] = v[1]; 
	
	string type = Util::str_to_lower(euler_type);
	if ( type != "2d") params["tz"] = v[2];
	
	float scale = get_scale();
	params["scale"] = scale;
	
	bool mirror = get_mirror();
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
	string euler_type;
	
	if (!rotation.has_key_ci("type") ){
			throw InvalidParameterException("argument dictionary does not contain the type entry");
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

	string type = Util::str_to_lower(euler_type);
	if (type == "2d") {
		validate_and_set_type(TWOD);
		az  = 0;
		alt = 0;
		phi = (float)rotation["alpha"] ;
	} else if ( type == "eman" ) {
		validate_and_set_type(THREED);
		az  = (float)rotation["az"] ;
		alt = (float)rotation["alt"]  ;
		phi = (float)rotation["phi"] ;
	} else if ( type == "imagic" ) {
		validate_and_set_type(THREED);
		az  = (float)rotation["alpha"] ;
		alt = (float)rotation["beta"]  ;
		phi = (float)rotation["gamma"] ;
	} else if ( type == "spider" ) {
		validate_and_set_type(THREED);
		az =  (float)rotation["phi"]    + 90.0f;
		alt = (float)rotation["theta"] ;
		phi = (float)rotation["psi"]    - 90.0f;
	} else if ( type == "xyz" ) {
		validate_and_set_type(THREED);
		cxtilt = cos( (M_PI/180.0f)*(float)rotation["xtilt"]);
		sxtilt = sin( (M_PI/180.0f)*(float)rotation["xtilt"]);
		cytilt = cos( (M_PI/180.0f)*(float)rotation["ytilt"]);
		sytilt = sin( (M_PI/180.0f)*(float)rotation["ytilt"]);	
		az =  (180.0f/M_PI)*atan2(-cytilt*sxtilt,sytilt)   + 90.0f ;
		alt = (180.0f/M_PI)*acos(cytilt*cxtilt)  ;
		phi = (float)rotation["ztilt"] +(180.0f/M_PI)*atan2(sxtilt,cxtilt*sytilt)   - 90.0f ;
	} else if ( type == "mrc" ) {
		validate_and_set_type(THREED);
		az  = (float)rotation["phi"]   + 90.0f ;
		alt = (float)rotation["theta"] ;
		phi = (float)rotation["omega"] - 90.0f ;
	} else if ( type == "quaternion" ) {
		validate_and_set_type(THREED);
		is_quaternion = 1;
		e0 = (float)rotation["e0"];
		e1 = (float)rotation["e1"];
		e2 = (float)rotation["e2"];
		e3 = (float)rotation["e3"];
	} else if ( type == "spin" ) {
		validate_and_set_type(THREED);
		is_quaternion = 1;
		Omega = (float)rotation["Omega"];
		e0 = cos(Omega*M_PI/360.0f);
		e1 = sin(Omega*M_PI/360.0f)* (float)rotation["n1"];
		e2 = sin(Omega*M_PI/360.0f)* (float)rotation["n2"];
		e3 = sin(Omega*M_PI/360.0f)* (float)rotation["n3"];
	} else if ( type == "sgirot" ) {
		validate_and_set_type(THREED);
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
		transform_type = UNKNOWN;
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

Dict Transform::get_rotation(const string& euler_type) const
{
	Dict result;
	
	float max = 1 - ERR_LIMIT;
	float scale;
	bool x_mirror;
	get_scale_and_mirror(scale,x_mirror);
	float cosalt=matrix[2][2]/scale;
	float x_mirror_scale = (x_mirror ? -1.0f : 1.0f);
	float inv_scale = 1.0f/scale;

	float az=0;
	float alt = 0;
	float phi=0;
	float phiS = 0; // like az   (but in SPIDER ZXZ)
	float psiS =0;  // like phi  (but in SPIDER ZYZ)

	// get alt, az, phi in EMAN convention
	if (cosalt > max) {  // that is, alt close to 0
		alt = 0;
		az=0;
		phi = (float)EMConsts::rad2deg*(float)atan2(x_mirror_scale*matrix[0][1], x_mirror_scale*matrix[0][0]); 
	}
	else if (cosalt < -max) { // alt close to pi
		alt = 180;
		az=0; 
		phi=360.0f-(float)EMConsts::rad2deg*(float)atan2(x_mirror_scale*matrix[0][1], x_mirror_scale*matrix[0][0]);
	}
	else {
		alt = (float)EMConsts::rad2deg*(float) acos(cosalt);
		az  = 360.0f+(float)EMConsts::rad2deg*(float)atan2(matrix[2][0], -matrix[2][1]);
		phi = 360.0f+(float)EMConsts::rad2deg*(float)atan2(x_mirror_scale*matrix[0][2], matrix[1][2]);
	}
// 	az=fmod(az+180.0f,360.0f)-180.0f;
// 	phi=fmod(phi+180.0f,360.0f)-180.0f;
	
	if (phi > 0) phi = fmod(phi,360.f);
	else phi = 360.f - fmod(fabs(phi),360.f);
	if (phi == 360.f) phi = 0.f;
	
	if (az > 0 ) az = fmod(az,360.f);
	else az = 360.f - fmod(fabs(az),360.f);
	if (az == 360.f) az = 0.f;
	
//   get phiS, psiS ; SPIDER
	if (fabs(cosalt) > max) {  // that is, alt close to 0
		phiS=0;
		psiS = phi;
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
			assert_consistent_type(TWOD);
			result["alpha"]  = phi;
	} else if (type == "eman") {
			assert_consistent_type(THREED);
			result["az"]  = az;
			result["alt"] = alt;
			result["phi"] = phi;
	} else if (type == "imagic") {
			assert_consistent_type(THREED);
			result["alpha"] = az;
			result["beta"] = alt;
			result["gamma"] = phi;
	} else if (type == "spider") {
		assert_consistent_type(THREED);
		result["phi"]   = phiS;  // The first Euler like az
		result["theta"] = alt;
		result["psi"]   = psiS;
	} else if (type == "mrc") {
		assert_consistent_type(THREED);
		result["phi"]   = phiS;
		result["theta"] = alt;
		result["omega"] = psiS;
	} else if (type == "xyz") {
		assert_consistent_type(THREED);
		xtilt = atan2(-sin((M_PI/180.0f)*phiS)*sin((M_PI/180.0f)*alt),cos((M_PI/180.0f)*alt));
		ytilt = asin(  cos((M_PI/180.0f)*phiS)*sin((M_PI/180.0f)*alt));
		ztilt = psiS*M_PI/180.0f - atan2(sin(xtilt), cos(xtilt) *sin(ytilt));

		xtilt=fmod(xtilt*180/M_PI+540.0f,360.0f) -180.0f;
		ztilt=fmod(ztilt*180/M_PI+540.0f,360.0f) -180.0f;

		result["xtilt"]  = xtilt;
		result["ytilt"]  = ytilt*180/M_PI;
		result["ztilt"]  = ztilt;
	} else if (type == "quaternion") {
		assert_consistent_type(THREED);
		result["e0"] = cosOover2 ;
		result["e1"] = sinOover2 * n1 ;
		result["e2"] = sinOover2 * n2;
		result["e3"] = sinOover2 * n3;
	} else if (type == "spin") {
		assert_consistent_type(THREED);
		result["Omega"] =360.0f* acos(cosOover2)/ M_PI ;
		result["n1"] = n1;
		result["n2"] = n2;
		result["n3"] = n3;
	} else if (type == "sgirot") {
		assert_consistent_type(THREED);
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
	if ( z != 0.0 ) {
		validate_and_set_type(THREED);
	}
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
	return v;
}

Vec2f Transform::get_trans_2d() const
{
// 	assert_consistent_type(TWOD);
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
	
	float c = ceilf(new_scale);
	float f = floorf(new_scale);
	float n_scale = new_scale;
	if (fabs(n_scale - c) < ERR_LIMIT) n_scale = c;
	else if (fabs(n_scale - f) < ERR_LIMIT) n_scale = f;
	
	
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
	
	
	float c = ceilf(scale);
	float f = floorf(scale);
	if (fabs(scale - c) < ERR_LIMIT) scale = c;
	else if (fabs(scale - f) < ERR_LIMIT) scale = f;
	
	return scale;
}

void Transform::orthogonalize()
{
	float inv_scale = 1.0f/get_scale();
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			matrix[i][j] *= inv_scale;
		}
	}
	set_mirror(false);
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
	
	float c = ceilf(scale);
	float f = floorf(scale);
	if (fabs(scale - c) < ERR_LIMIT) scale = c;
	else if (fabs(scale - f) < ERR_LIMIT) scale = f;
}

float Transform::get_determinant() const
{
	float det = matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[2][1]*matrix[1][2]);
	det -= matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[2][0]*matrix[1][2]);
	det += matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[2][0]*matrix[1][1]);

	float c = ceilf(det);
	float f = floorf(det);
	if (fabs(det - c) < ERR_LIMIT) det = c;
	else if (fabs(det - f) < ERR_LIMIT) det = f;
	
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

Transform Transform::transpose() const
{
	Transform t(*this);
	t.transpose_inplace();
	return t;
}

void Transform::transpose_inplace() {
	float tempij;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < i; j++) {
			tempij= matrix[i][j];
			matrix[i][j] = matrix[j][i];
			matrix[j][i] = tempij;
		}
	}
}

void Transform::validate_and_set_type(const Transform::TransformType type ) {
	if (transform_type == UNKNOWN ) {
		transform_type = type; // This is fine this is the first time the type is being set
	} else if (transform_type != type ) {
		string message = "can't treat a " + transform_type_to_string(transform_type) + " Transform as though it were a " +  transform_type_to_string(type) + " Transform type";
		throw UnexpectedBehaviorException("Error, " + message);
	} 
}

void Transform::assert_consistent_type(const Transform::TransformType type )  const{
	if (transform_type != UNKNOWN ) {// if the type is unknown then maybe it's the identity and the user is asking for parameters?
		if (transform_type != type ) {
			string message = "can't treat a " + transform_type_to_string(transform_type) + " Transform type as though it were a " +  transform_type_to_string(type) + " Transform type";
			throw UnexpectedBehaviorException("Error, " + message);
		}
	}
}

string Transform::transform_type_to_string(const Transform::TransformType type ) const {
	if ( type == UNKNOWN ) {
		return string("UNKNOWN");
	} else if ( type == TWOD ) {
		return string("2D");
	}
	else if ( type == THREED ) {
		return string("3D");
	} else {
		throw InvalidParameterException("Error, unknown type in Transform::transform_type_to_string");
	}
}

Transform EMAN::operator*(const Transform & M2, const Transform & M1)     // YYY
{
	if (M2.transform_type != Transform::UNKNOWN && M1.transform_type != Transform::UNKNOWN ) {
		if (M2.transform_type != M1.transform_type) {
			string message = "can't multiply a " + M2.transform_type_to_string(M2.transform_type) + " Transform type against a " + M1.transform_type_to_string(M1.transform_type) + " Transform type";
			throw UnexpectedBehaviorException("Error, " + message);
		}
	}
	Transform result;
	for (int i=0; i<3; i++) {
		for (int j=0; j<4; j++) {
			result[i][j] = M2[i][0] * M1[0][j] +  M2[i][1] * M1[1][j] + M2[i][2] * M1[2][j];
		}
		result[i][3] += M2[i][3];
	}
	
	return result;
}

void Transform::assert_valid_2d() {
	
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
		psiS = phi;
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



template <> Factory < Symmetry3D >::Factory()
{
	force_add(&CSym::NEW);
	force_add(&DSym::NEW);
	force_add(&HSym::NEW);
	force_add(&TetrahedralSym::NEW);
	force_add(&OctahedralSym::NEW);
	force_add(&IcosahedralSym::NEW);
}



void EMAN::dump_symmetries()
{
	dump_factory < Symmetry3D > ();
}

map<string, vector<string> > EMAN::dump_symmetries_list()
{
	return dump_factory_list < Symmetry3D > ();
}

template <>
Symmetry3D* Factory < Symmetry3D >::get(const string & instancename)
{
	init();
	
	unsigned int n = instancename.size();
	if ( n == 0 ) throw NotExistingObjectException(instancename, "Empty instance name!");
	
	char leadingchar = instancename[0];
	if (leadingchar == 'c' || leadingchar == 'd' || leadingchar == 'h' ) {
		Dict parms;
		if (n > 1) {
			int nsym = atoi(instancename.c_str() + 1);
			parms["nsym"] = nsym;
		}
		
		if (leadingchar == 'c') {
			return get("c",parms);	
		}
		if (leadingchar == 'd') {
			return get("d",parms);	
		}
		if (leadingchar == 'h') {
			return get("h",parms);	
		}
		
// 		delete lc;
	}
	else if ( instancename == "icos" || instancename == "oct" || instancename == "tet" )
	{
		map < string, InstanceType >::iterator fi =
				my_instance->my_dict.find(instancename);
		if (fi != my_instance->my_dict.end()) {
			return my_instance->my_dict[instancename] ();
		}
	
		string lower = instancename;
		for (unsigned int i=0; i<lower.length(); i++) lower[i]=tolower(lower[i]);
	
		fi = my_instance->my_dict.find(lower);
		if (fi != my_instance->my_dict.end()) {
			return my_instance->my_dict[lower] ();
		}
		
		throw NotExistingObjectException(instancename, "No such an instance existing");
	}
	else throw NotExistingObjectException(instancename, "No such an instance existing");
	
	throw NotExistingObjectException(instancename, "No such an instance existing");
}

template <> Factory < OrientationGenerator >::Factory()
{
	force_add(&EmanOrientationGenerator::NEW);
	force_add(&RandomOrientationGenerator::NEW);
	force_add(&EvenOrientationGenerator::NEW);
	force_add(&SaffOrientationGenerator::NEW);
	force_add(&OptimumOrientationGenerator::NEW);
}



void EMAN::dump_orientgens()
{
	dump_factory < OrientationGenerator > ();
}

map<string, vector<string> > EMAN::dump_orientgens_list()
{
	return dump_factory_list < OrientationGenerator > ();
}

vector<Transform> Symmetry3D::gen_orientations(const string& generatorname, const Dict& parms)
{
	ENTERFUNC;
	vector<Transform> ret;
	OrientationGenerator *g = Factory < OrientationGenerator >::get(generatorname, parms);
	if (g) {
		ret = g->gen_orientations(this);
		if( g )
		{
			delete g;
			g = 0;
		}
	}
	else throw;
	
	EXITFUNC;
	
	return ret;
}

int EmanOrientationGenerator::get_orientations_tally(const Symmetry3D* const sym, const float& delta) const
{
	//FIXME THIS IS SO SIMILAR TO THE gen_orientations function that they should be probably use
	// a common routine - SAME ISSUE FOR OTHER ORIENTATION GENERATORS
	bool inc_mirror = params.set_default("inc_mirror",false);
	Dict delimiters = sym->get_delimiters(inc_mirror);
	float altmax = delimiters["alt_max"];
	float azmax = delimiters["az_max"];

	float alt_iterator = 0.0f;
	
	// #If it's a h symmetry then the alt iterator starts at very close
	// #to the altmax... the object is a h symmetry then it knows its alt_min...
	if (sym->is_h_sym()) alt_iterator = delimiters["alt_min"];
	
	int tally = 0;
	while ( alt_iterator <= altmax ) {
		float h = get_az_delta(delta,alt_iterator, sym->get_max_csym() );

		// not sure what this does code taken from EMAN1 - FIXME original author add comments
		if ( (alt_iterator > 0) && ( (azmax/h) < 2.8) ) h = azmax / 2.1f;
		else if (alt_iterator == 0) h = azmax;
			
		float az_iterator = 0.0;
		
		float azmax_adjusted = azmax;
		
		// if this is odd c symmetry, and we're at the equator, and we're excluding the mirror then
		// half the equator is redundant (it is the mirror of the other half)
		if (sym->is_c_sym() && !inc_mirror && alt_iterator == altmax && (sym->get_nsym() % 2 == 1 ) ){
			azmax_adjusted /= 2.0;
		}
		// at the azimuthal boundary in c symmetry and tetrahedral symmetry we have come
		// full circle, we must not include it
		else if (sym->is_c_sym() || sym->is_tet_sym() ) {
			azmax_adjusted -=  h/4.0f;
		}
		// If we're including the mirror then in d and icos and oct symmetry the azimuthal
		// boundary represents coming full circle, so must be careful to exclude it
		else if (inc_mirror && ( sym->is_d_sym() || sym->is_platonic_sym() ) )  {
			azmax_adjusted -=  h/4.0f;
		}
		// else do nothing - this means that we're including the great arc traversing
		// the full range of permissable altitude angles at azmax.
		// This happens in d symmetry, and in the icos and oct symmetries, when the mirror
		// portion of the asymmetric unit is being excluded
		

		while ( az_iterator <= azmax_adjusted ) {
			// FIXME: add an intelligent comment - this was copied from old code	
			if ( az_iterator > 180.0 && alt_iterator > 180.0/(2.0-0.001) && alt_iterator < 180.0/(2.0+0.001) ) {
				az_iterator +=  h;
				continue;
			}
			
			if (sym->is_platonic_sym()) {
				if ( sym->is_in_asym_unit(alt_iterator, az_iterator,inc_mirror) == false ) {
					az_iterator += h;
					continue;
				}
			}
				
			tally++;
			if ( sym->is_h_sym() && inc_mirror && alt_iterator != (float) delimiters["alt_min"] ) {
				tally++;
			}
			az_iterator += h;
		}
		alt_iterator += delta;
	}
	return tally;
}

float OrientationGenerator::get_optimal_delta(const Symmetry3D* const sym, const int& n) const
{
	
	float delta_soln = 360.0f/sym->get_max_csym();
	float delta_upper_bound = delta_soln;
	float delta_lower_bound = 0.0;
	
	int prev_tally = -1;
	// This is an example of a divide and conquer approach, the possible values of delta are searched
	// like a binary tree
	
	bool soln_found = false;
	
	while ( soln_found == false ) {
		int tally = get_orientations_tally(sym,delta_soln);
		if ( tally == n ) soln_found = true;
		else if ( (delta_upper_bound - delta_lower_bound) < 0.0001 ) {
			// If this is the case, the requested number of projections is practically infeasible
			// in which case just return the nearest guess
			soln_found = true;
			delta_soln = (delta_upper_bound+delta_lower_bound)/2.0f;
		}
		else if (tally < n) {
			delta_upper_bound = delta_soln;
			delta_soln = delta_soln - (delta_soln-delta_lower_bound)/2.0f;
		}
		else  /* tally > n*/{
			delta_lower_bound = delta_soln;
			delta_soln = delta_soln  + (delta_upper_bound-delta_soln)/2.0f;
		}
		prev_tally = tally;
	}
	
	return delta_soln;
}

bool OrientationGenerator::add_orientation(vector<Transform>& v, const float& az, const float& alt) const
{
	bool randphi = params.set_default("random_phi",false);
	float phi = 0.0;
	if (randphi) phi = Util::get_frand(0.0f,359.99999f);
	float phitoo = params.set_default("phitoo",0.0f);
	if ( phitoo < 0 ) throw InvalidValueException(phitoo, "Error, if you specify phitoo is must be positive");
	Dict d;
	d["type"] = "eman";
	d["az"] = az;
	d["alt"] = alt;
	d["phi"] = phi;
	Transform t(d);
	v.push_back(t);
	if ( phitoo != 0 ) {
		if (phitoo < 0) return false;
		else {
			for ( float p = phitoo; p <= 360.0-phitoo; p+= phitoo )
			{
				d["phi"] = fmod(phi+p,360);
				Transform t(d);
				v.push_back(t);
			}
		}
	}
	return true;
}

float EmanOrientationGenerator::get_az_delta(const float& delta,const float& altitude, const int maxcsym) const
{
	// convert altitude into radians
	float tmp = (float)(EMConsts::deg2rad * altitude);
			
	// This is taken from EMAN1 project3d.C
	float h=floor(360.0f/(delta*1.1547f));	// the 1.1547 makes the overall distribution more like a hexagonal mesh
	h=(int)floor(h*sin(tmp)+.5f);
	if (h==0) h=1;
	h=abs(maxcsym)*floor(h/(float)abs(maxcsym)+.5f);
	if ( h == 0 ) h = (float)maxcsym;
	h=2.0f*M_PI/h;
	
	return (float)(EMConsts::rad2deg*h);
}

vector<Transform> EmanOrientationGenerator::gen_orientations(const Symmetry3D* const sym) const
{	
	float delta = params.set_default("delta", 0.0f);
	int n = params.set_default("n", 0);
	
	if ( delta <= 0 && n <= 0 ) throw InvalidParameterException("Error, you must specify a positive non-zero delta or n");
	if ( delta > 0 && n > 0 ) throw InvalidParameterException("Error, the delta and the n arguments are mutually exclusive");
	
	if ( n > 0 ) {
		delta = get_optimal_delta(sym,n);
	}
	
	bool inc_mirror = params.set_default("inc_mirror",false);
	Dict delimiters = sym->get_delimiters(inc_mirror);
	float altmax = delimiters["alt_max"];
	float azmax = delimiters["az_max"];
	
	bool perturb = params.set_default("perturb",false);
	
	float alt_iterator = 0.0;
	
	// #If it's a h symmetry then the alt iterator starts at very close
	// #to the altmax... the object is a h symmetry then it knows its alt_min...
	if (sym->is_h_sym()) alt_iterator = delimiters["alt_min"];
	
	vector<Transform> ret;
	while ( alt_iterator <= altmax ) {
		float h = get_az_delta(delta,alt_iterator, sym->get_max_csym() );

		// not sure what this does code taken from EMAN1 - FIXME original author add comments
		if ( (alt_iterator > 0) && ( (azmax/h) < 2.8) ) h = azmax / 2.1f;
		else if (alt_iterator == 0) h = azmax;
			
		float az_iterator = 0.0;
		
		float azmax_adjusted = azmax;
		
		// if this is odd c symmetry, and we're at the equator, and we're excluding the mirror then
		// half the equator is redundant (it is the mirror of the other half)
		if (sym->is_c_sym() && !inc_mirror && alt_iterator == altmax && (sym->get_nsym() % 2 == 1 ) ){
			azmax_adjusted /= 2.0;
		}
		// at the azimuthal boundary in c symmetry and tetrahedral symmetry we have come
		// full circle, we must not include it
		else if (sym->is_c_sym() || sym->is_tet_sym() ) {
			azmax_adjusted -=  h/4.0f;
		}
		// If we're including the mirror then in d and icos and oct symmetry the azimuthal
		// boundary represents coming full circle, so must be careful to exclude it
		else if (inc_mirror && ( sym->is_d_sym() || sym->is_platonic_sym() ) )  {
			azmax_adjusted -=  h/4.0f;
		}
		// else do nothing - this means that we're including the great arc traversing
		// the full range of permissable altitude angles at azmax.
		// This happens in d symmetry, and in the icos and oct symmetries, when the mirror
		// portion of the asymmetric unit is being excluded
		

		while ( az_iterator <= azmax_adjusted ) {
			// FIXME: add an intelligent comment - this was copied from old code	
// 			if ( az_iterator > 180.0 && alt_iterator > 180.0/(2.0-0.001) && alt_iterator < 180.0/(2.0+0.001) ) {
// 				az_iterator +=  h;
// 				continue;
// 			}
// 			// Now that I am handling the boundaries very specifically, I don't think we need
			// the above if statement. But I am leaving it there in case I need to reconsider.
			
			float alt_soln = alt_iterator;
			float az_soln = az_iterator;
			
			if (sym->is_platonic_sym()) {
				if ( sym->is_in_asym_unit(alt_soln, az_soln,inc_mirror) == false ) {
					az_iterator += h;
					continue;
				}
				
				// Some objects have alignment offsets (icos and tet particularly)
				az_soln += sym->get_az_alignment_offset();
			}
				
			if ( perturb &&  alt_soln != 0 ) {
				alt_soln += Util::get_gauss_rand(0.0f,.25f*delta);
				az_soln += Util::get_gauss_rand(0.0f,h/4.0f);
			}
				
			add_orientation(ret,az_soln,alt_soln);
		
			// Add helical symmetry orientations on the other side of the equator (if we're including
			// mirror orientations) 
			if ( sym->is_h_sym() && inc_mirror && alt_iterator != (float) delimiters["alt_min"] ) {
				add_orientation(ret, az_soln,2.0f*(float)delimiters["alt_min"]-alt_soln);
			}
			az_iterator += h;
			
		}
		alt_iterator += delta;
	}
	
	return ret;
}

vector<Transform> RandomOrientationGenerator::gen_orientations(const Symmetry3D* const sym) const
{
	int n = params.set_default("n", 0);
	
	if ( n <= 0 ) throw InvalidParameterException("You must specify a positive, non zero n for the Random Orientation Generator");
	
	bool phitoo = params.set_default("phitoo", false);
	bool inc_mirror = params.set_default("inc_mirror", false);
	
	vector<Transform> ret;
	
	int i = 0;
	Dict d("type","eman");
	while ( i < n ){
		float u1 =  Util::get_frand(-1.0f,1.0f);
		float u2 =  Util::get_frand(-1.0f,1.0f);
		float s = u1*u1 + u2*u2;
		if ( s > 1.0 ) continue;
		float alpha = 2.0f*sqrtf(1.0f-s);
		float x = alpha * u1;
		float y = alpha * u2;
		float z = 2.0f*s-1.0f;
		
		float altitude = (float)EMConsts::rad2deg*acos(z);
		float azimuth = (float)EMConsts::rad2deg*atan2(y,x);
		
		float phi = 0.0;
		if ( phitoo ) phi = Util::get_frand(0.0,359.9999);
		
		d["az"] = azimuth; d["phi"] = phi; d["alt"] = altitude;
		Transform t(d);
		
		if ( !(sym->is_c_sym() && sym->get_nsym() == 1)) t = sym->reduce(t); //reduce doesn't make sense for C1 symmetry
		
		if ( !sym->is_in_asym_unit(altitude,azimuth,inc_mirror) ){
			// is_in_asym_unit has returned the wrong value!
			// FIXME
// 			cout << "warning, there is an unresolved issue - email D Woolford" << endl;
		}
		ret.push_back(t);
		i++;
	}
	return ret;
}

int EvenOrientationGenerator::get_orientations_tally(const Symmetry3D* const sym, const float& delta) const
{
	bool inc_mirror = params.set_default("inc_mirror",false);
	Dict delimiters = sym->get_delimiters(inc_mirror);
	float altmax = delimiters["alt_max"];
	float azmax = delimiters["az_max"];
	
	float altmin = 0.0;
	// #If it's a h symmetry then the alt iterator starts at very close
	// #to the altmax... the object is a h symmetry then it knows its alt_min...
	if (sym->is_h_sym()) altmin = delimiters["alt_min"];
	
	int tally = 0;
	
	for (float alt = altmin; alt <= altmax; alt += delta) {
		float detaz;
		int lt;
		if ((0.0 == alt)||(180.0 == alt)) {
			detaz = 360.0f;
			lt = 1;
		} else {
			detaz = delta/(float)sin(alt*EMConsts::deg2rad);
			lt = int(azmax/detaz)-1;
			if (lt < 1) lt = 1;
			detaz = azmax/(float)lt;
		}
		for (int i = 0; i < lt; i++) {
			float az = (float)i*detaz;
			if (sym->is_platonic_sym()) {
				if ( sym->is_in_asym_unit(alt, az,inc_mirror) == false ) continue;
			}
			tally++;
			if ( sym->is_h_sym() && inc_mirror && alt != altmin ) {
				tally++;
			}
		}
	}

	return tally;	
}

vector<Transform> EvenOrientationGenerator::gen_orientations(const Symmetry3D* const sym) const
{
	float delta = params.set_default("delta", 0.0f);
	int n = params.set_default("n", 0);
	
	if ( delta <= 0 && n <= 0 ) throw InvalidParameterException("Error, you must specify a positive non-zero delta or n");
	if ( delta > 0 && n > 0 ) throw InvalidParameterException("Error, the delta and the n arguments are mutually exzclusive");
	
	if ( n > 0 ) {
		delta = get_optimal_delta(sym,n);
	}
	
	bool inc_mirror = params.set_default("inc_mirror",false);
	Dict delimiters = sym->get_delimiters(inc_mirror);
	float altmax = delimiters["alt_max"];
	float azmax = delimiters["az_max"];
	
	float altmin = 0.0;
	// If it's a h symmetry then the alt iterator starts at very close
	// to the altmax... the object is a h symmetry then it knows its alt_min...
	if (sym->is_h_sym()) altmin = delimiters["alt_min"];
	
	vector<Transform> ret;
	
	for (float alt = altmin; alt <= altmax; alt += delta) {
		float detaz;
		int lt;
		if ((0.0 == alt)||(180.0 == alt)) {
			detaz = 360.0f;
			lt = 1;
		} else {
			detaz = delta/(float)sin(alt*EMConsts::deg2rad);
			lt = int(azmax/detaz)-1;
			if (lt < 1) lt = 1;
			detaz = azmax/(float)lt;
		}
		for (int i = 0; i < lt; i++) {
			float az = (float)i*detaz;
			if (sym->is_platonic_sym()) {
				if ( sym->is_in_asym_unit(alt, az,inc_mirror) == false ) continue;
			}
			add_orientation(ret,az,alt);
			if ( sym->is_h_sym() && inc_mirror && alt != altmin ) {
				add_orientation(ret,az,2.0f*altmin-alt);
			}
		}
	}

	return ret;
}

int SaffOrientationGenerator::get_orientations_tally(const Symmetry3D* const sym, const float& delta) const
{
	bool inc_mirror = params.set_default("inc_mirror",false);
	Dict delimiters = sym->get_delimiters(inc_mirror);
	float altmax = delimiters["alt_max"];
	float azmax = delimiters["az_max"];
	
	float altmin = 0.0;
	// #If it's a h symmetry then the alt iterator starts at very close
	// #to the altmax... the object is a h symmetry then it knows its alt_min...
	if (sym->is_h_sym()){
		altmin = delimiters["alt_min"];
		if (inc_mirror) {
			altmin -= (float) sym->get_params()["equator_range"];
		}
	}
	
	float Deltaz = (float)(cos(altmax*EMConsts::deg2rad)-cos(altmin*EMConsts::deg2rad));
	float s = delta*M_PI/180.0f;
	float NFactor = 3.6f/s;
	float wedgeFactor = fabs( Deltaz*(azmax)/720.0f) ;
	int NumPoints   =  static_cast<int> (NFactor*NFactor*wedgeFactor);
	
	int tally = 0;
	if (!sym->is_h_sym()) ++tally;
	float az = 0.0f;
	float dz = (float)cos(altmin*EMConsts::deg2rad);
	for(int i = 1; i < NumPoints; ++i ){
		float z = dz + Deltaz* (float)i/ float(NumPoints-1);
		float r= sqrt(1.0f-z*z);
		az = fmod(az + delta/r,azmax);
		float alt = (float)(acos(z)*EMConsts::rad2deg);
		if (sym->is_platonic_sym()) {
			if ( sym->is_in_asym_unit(alt,az,inc_mirror) == false ) continue;
		}
		tally++;
	}

	return tally;
}

vector<Transform> SaffOrientationGenerator::gen_orientations(const Symmetry3D* const sym) const
{
	float delta = params.set_default("delta", 0.0f);
	int n = params.set_default("n", 0);
	
	if ( delta <= 0 && n <= 0 ) throw InvalidParameterException("Error, you must specify a positive non-zero delta or n");
	if ( delta > 0 && n > 0 ) throw InvalidParameterException("Error, the delta and the n arguments are mutually exclusive");
	
	if ( n > 0 ) {
		delta = get_optimal_delta(sym,n);
	}
	
// 	if ( sym->is_platonic_sym() ) return gen_platonic_orientations(sym, delta);
	
	bool inc_mirror = params.set_default("inc_mirror",false);
	Dict delimiters = sym->get_delimiters(inc_mirror);
	float altmax = delimiters["alt_max"];
	float azmax = delimiters["az_max"];
	
	float altmin = 0.0;
	// #If it's a h symmetry then the alt iterator starts at very close
	// #to the altmax... the object is a h symmetry then it knows its alt_min...
	if (sym->is_h_sym()){
		 altmin = delimiters["alt_min"];
		if (inc_mirror) {
			altmin -= (float) sym->get_params()["equator_range"];
		}
	}
	
	float Deltaz = (float)(cos(altmax*EMConsts::deg2rad)-cos(altmin*EMConsts::deg2rad));
	float s = delta*M_PI/180.0f;
	float NFactor = 3.6f/s;
	float wedgeFactor = fabs( Deltaz*(azmax)/720.0f) ;
	int NumPoints   =  static_cast<int> (NFactor*NFactor*wedgeFactor);
	
	vector<Transform> ret;
	
	if (!sym->is_h_sym()) add_orientation(ret,0,0);
	float az = 0.0;
	float dz = (float)cos(altmin*EMConsts::deg2rad);
	for(int i = 1; i < NumPoints; ++i ){
		float z = dz + Deltaz* (float)i/ float(NumPoints-1);
		float r= sqrt(1.0f-z*z);
		az = fmod(az + delta/r,azmax);
		float alt = (float)(acos(z)*EMConsts::rad2deg);
		if (sym->is_platonic_sym()) {
			if ( sym->is_in_asym_unit(alt,az,inc_mirror) == false ) continue;
		}
		add_orientation(ret,az,alt);
	}

	return ret;
}

int OptimumOrientationGenerator::get_orientations_tally(const Symmetry3D* const sym, const float& delta) const
{
	string deltaoptname = params.set_default("use","saff");
	Dict a;
	a["inc_mirror"] = (bool)params.set_default("inc_mirror",false);
	OrientationGenerator *g = Factory < OrientationGenerator >::get(deltaoptname,a);
	if (g) {
		int tally = g->get_orientations_tally(sym,delta);
		delete g;
		g = 0;
		return tally;
	}
	else throw;
}

vector<Transform> OptimumOrientationGenerator::gen_orientations(const Symmetry3D* const sym) const
{	
	float delta = params.set_default("delta", 0.0f);
	int n = params.set_default("n", 0);
	
	bool inc_mirror = params.set_default("inc_mirror",false);
	
	if ( delta <= 0 && n <= 0 ) throw InvalidParameterException("Error, you must specify a positive non-zero delta or n");
	if ( delta > 0 && n > 0 ) throw InvalidParameterException("Error, the delta and the n arguments are mutually exclusive");
	
	string generatorname = params.set_default("use","saff");
	
	if ( n > 0 && generatorname != RandomOrientationGenerator::NAME ) {
		params["delta"] = get_optimal_delta(sym,n);
		params["n"] = (int)0;
	}
	
	// Force the orientation generator to include the mirror - this is because 
	// We will enventually use it to generate orientations over the intire sphere
	// which is C1 symmetry, with the inc_mirror flag set to true
	params["inc_mirror"] = true;
	OrientationGenerator* g = Factory < OrientationGenerator >::get(generatorname);
	g->set_params(copy_relevant_params(g));

	
	// get the starting orientation distribution
	CSym* unit_sphere = new CSym();
	Dict nsym; nsym["nsym"] = 1; unit_sphere->set_params(nsym);
	
	vector<Transform> unitsphereorientations = g->gen_orientations(unit_sphere);
	delete g; g = 0;
	delete unit_sphere; unit_sphere = 0;
	
	vector<Vec3f> angles = optimize_distances(unitsphereorientations);
	
	vector<Transform> ret;
	for (vector<Vec3f>::const_iterator it = angles.begin(); it != angles.end(); ++it ) {
		if ( sym->is_in_asym_unit((*it)[1],(*it)[0],inc_mirror) ) {
			add_orientation(ret,(*it)[0],(*it)[1]);
		}
	}
	
	// reset the params to what they were before they were acted upon by this class
	params["inc_mirror"] = inc_mirror;
	params["delta"] = delta;
	params["n"] = n;
	
	return ret;
}

vector<Vec3f> OptimumOrientationGenerator::optimize_distances(const vector<Transform>& v) const
{
	vector<Vec3f> points;
	
	for (vector<Transform>::const_iterator it = v.begin(); it != v.end(); ++it ) {
		points.push_back(Vec3f(0,0,1)*(*it));
	}
	
	if ( points.size() >= 2 ) {
		int max_it = 1000;
		float percentage = 0.01f;
		
		for ( int i = 0; i < max_it; ++i ){
			unsigned int p1 = 0;
			unsigned int p2 = 1;
			
			float distsquared = (points[p1]-points[p2]).squared_length();
			
			// Find the nearest points
			for(unsigned int j = 0; j < points.size(); ++j) {
				for(unsigned int k = j+1; k < points.size(); ++k) {
					float d = (points[j]-points[k]).squared_length();
					if ( d < distsquared ) {
						distsquared = d;
						p1 = j;
						p2 = k;
					}
				}
			}
			
			// Move them apart by a small fraction
			Vec3f delta = percentage*(points[p2]-points[p1]);
			
			points[p2] += delta;
			points[p2].normalize();
			points[p1] -= delta;
			points[p1].normalize();
		}
	}
	
	vector<Vec3f> ret;
	for (vector<Vec3f>::const_iterator it = points.begin(); it != points.end(); ++it ) {
		float altitude = (float)(EMConsts::rad2deg*acos((*it)[2]));
		float azimuth = (float)(EMConsts::rad2deg*atan2((*it)[1],(*it)[0]));
		ret.push_back(Vec3f(90.0f+azimuth,altitude,0));
	}
	
	return ret;
}
// THIS IS DWOOLFORDS FIRST SHOT AT EXTRACTING PHIL'S PLATONIC STUFF FROM SPARX
// It didn't work.
// vector<Transform3D> SaffOrientationGenerator::gen_platonic_orientations(const Symmetry3D* const sym, const float& delta) const	
// {
// 	float scrunch = 0.9; //closeness factor to eliminate oversampling corners
// 		
// 	float fudge; //# fudge is a factor used to adjust phi steps
// 	float m = static_cast<float>(sym->get_max_csym());
// 	if ( sym->get_name() == TetrahedralSym::NAME ) fudge=0.9;
// 	else if ( sym->get_name() == OctahedralSym::NAME ) fudge=0.8;
// 	else if ( sym->get_name() == IcosahedralSym::NAME) fudge=0.95;
// 	else throw; // this should not happen
// 		
// 	float n=3.0;
// 	float OmegaR = 2.0*M_PI/m;
// 	float cosOmega= cos(OmegaR);
// 	int Edges  = static_cast<int>(2.0*m*n/(2.0*(m+n)-m*n));
// 	int Faces  = static_cast<int>(2*Edges/n);
// 	float Area   = 4*M_PI/Faces/3.0; // also equals  2*pi/3 + Omega
// 	float costhetac = cosOmega/(1-cosOmega);
// 	float deltaRad= delta*M_PI/180.0;
// 	
// 	int	NumPoints = static_cast<int>(Area/(deltaRad*deltaRad));
// 	float fheight = 1.0f/sqrt(3.0f)/(tan(OmegaR/2.0f));
// 	float z0      = costhetac; // initialize loop
// 	float z       = z0;
// 	float phi     = 0;
// 	float Deltaz  = (1-costhetac);
// 	
// 	vector<Transform3D> ret;
// 	ret.push_back(Transform3D(phi, acos(z)*EMConsts::rad2deg, 0 ));
// 
// 	vector<Vec3f> points;
// 	points.push_back(Vec3f(sin(acos(z))*cos(phi*EMConsts::deg2rad) ,  sin(acos(z))*sin(phi*EMConsts::deg2rad) , z) );
// 	//nLast=  [ sin(acos(z))*cos(phi*piOver) ,  sin(acos(z))*sin(phi*piOver) , z]
// 	//nVec.append(nLast)	
// 		
// 	for(int k = 0; k < NumPoints-1; ++k ) {
// 		z = z0 + Deltaz*(float)k/(float)(NumPoints-1);
// 		float r= sqrt(1-z*z);
// 		float phiMax =180.0*OmegaR/M_PI/2.0;
// 		// Is it higher than fhat or lower
// 		if (z<= fheight && false) {
// 			float thetaR   = acos(z); 
// 			float cosStuff = (cos(thetaR)/sin(thetaR))*sqrt(1. - 2 *cosOmega);
// 			phiMax   =  180.0*( OmegaR - acos(cosStuff))/M_PI;
// 		}
// 		float angleJump = fudge* delta/r;
// 		phi = fmod(phi + angleJump,phiMax);
// // 		anglesNew = [phi,180.0*acos(z)/pi,0.];
// // 		Vec3f nNew( sin(acos(z))*cos(phi*EMConsts::deg2rad) ,  sin(acos(z))*sin(phi*EMConsts::deg2rad) , z);
// // 		float mindiff = acos(nNew.dot(points[0]));
// // 		for(unsigned int l = 0; l < points.size(); ++l ) {
// // 			float dx = acos(nNew.dot(points[l]));
// // 			if (dx < mindiff ) mindiff = dx;
// // 		}
// // 		if (mindiff > (angleJump*EMConsts::deg2rad *scrunch) ){
// // 			points.push_back(nNew);
// 			cout << "phi " << phi << " alt " << acos(z)*EMConsts::rad2deg << " " << z << endl;
// 			ret.push_back(Transform3D(phi, acos(z)*EMConsts::rad2deg, 0 ));
// // 		}
// 	}
// 	ret.push_back(Transform3D(0,0,0 ));
// 	
// 	return ret;
// }



void verify(const Vec3f& tmp, float * plane, const string& message )
{
	cout << message << " residual " << plane[0]*tmp[0]+plane[1]*tmp[1]+plane[2]*tmp[2] + plane[3]  << endl;
}

Transform Symmetry3D::reduce(const Transform& t, int n) const
{
	// Determine which asym unit the given asym unit is in
	int soln = in_which_asym_unit(t);
	
	// This should never happen
	if ( soln == -1 ) {
		cout << "error, no solution found!" << endl;
		throw;
	}
	
	// Get the symmetry operation corresponding to the intersection asymmetric unit
	Transform nt = get_sym(soln);
	// Transpose it (invert it)
	nt.transpose_inplace();
	// Now we can transform the argument orientation into the default asymmetric unit
	nt  = t*nt;
	// Now that we're at the default asymmetric unit, we can map into the requested asymmunit by doing this
	if ( n != 0 ) {
		nt = nt*get_sym(n);
	}
	// Done!
	return nt;
	
}

int Symmetry3D::in_which_asym_unit(const Transform& t3d) const
{
	// Here it is assumed that final destination of the orientation (as encapsulated in the t3d object) is
	// in the z direction, so in essence we will start in the direction z and 'undo' the orientation to get the real
	// direction
	Vec3f p(0,0,1);
	
	Transform o(t3d);
	// Orientations are alway transposed when dealing with asymmetric units, projections,etc
	// We take the transpose to 'undo' the transform and get the true direction of the point.
	o.transpose_inplace();
	// Figure out where the point would end up. No we could just as easily not transpose and do
	// left multiplation (as in what occurs in the FourierReconstructor during slice insertion)
	p = o*p;
	
	// This is empty space for storing the equation of a plane
	float* plane = new float[4];
	
	// Get a set of triangles that, if projected through, covers the entire asymmetric unit.
	// Sometimes this will only be one triangle, sometimes it's two, occassionally it's 4.
	vector<vector<Vec3f> >triangles = get_asym_unit_triangles(true);
	typedef vector<vector<Vec3f> >::const_iterator cit;
	
	for(int i = 0; i < get_nsym(); ++i) {
		
		for( cit it = triangles.begin(); it != triangles.end(); ++it )
		{
			// For each given triangle
			vector<Vec3f> points = *it;
			if ( i != 0 ) {
				for (vector<Vec3f>::iterator it = points.begin(); it != points.end(); ++it ) {
					// Rotate the points in the triangle so that the triangle occupies the
					// space of the current asymmetric unit
					*it = (*it)*get_sym(i); 
				}
			}
			
			// Determine the equation of the plane for the points, store it in plane
			Util::equation_of_plane(points[0],points[2],points[1],plane);

			Vec3f tmp = p;
	
			// Determine the intersection of p with the plane - do this by finding out how much p should be scaled by
			float scale = plane[0]*tmp[0]+plane[1]*tmp[1]+plane[2]*tmp[2];
			if ( scale != 0 )
				scale = -plane[3]/scale;
			else {
				// parralel!
				continue;
			}
			
			// If the scale factor is less than zero, then p is definitely not in this asymmetric unit 
			if (scale <= 0) continue;	
			
			// This is the intersection point
			Vec3f pp = tmp*scale;
			
			// Now we have to see if the point p is inside the region bounded by the points, or if it is outside
			// If it is inside the region then p is in this asymmetric unit.
			
			// This formula take from FIXME fill in once I get to work
			Vec3f v = points[2]-points[0];
			Vec3f u = points[1]-points[0];
			Vec3f w = pp - points[0];
			
			float udotu = u.dot(u);
			float udotv = u.dot(v); 
			float udotw = u.dot(w);
			float vdotv = v.dot(v);
			float vdotw = v.dot(w);
			
			float d = 1.0f/(udotv*udotv - udotu*vdotv);
			float s = udotv*vdotw - vdotv*udotw;
			s *= d;
			
			float t = udotv*udotw - udotu*vdotw;
			t *= d;
			
			// We've done a few multiplications, so detect when there are tiny residuals that may throw off the final 
			// decision
			if (fabs(s) < Transform::ERR_LIMIT ) s = 0;
			if (fabs(t) < Transform::ERR_LIMIT ) t = 0;
			
			if ( fabs((fabs(s)-1.0)) < Transform::ERR_LIMIT ) s = 1;
			if ( fabs((fabs(t)-1.0)) < Transform::ERR_LIMIT ) t = 1;
			
			// The final decision, if this is true then we've hit the jackpot
			if ( s >= 0 && t >= 0 && (s+t) <= 1 ) {
				delete [] plane;
				return i;
			}
		}
	}

	delete [] plane;
	return -1;
}

vector<Transform> Symmetry3D::get_touching_au_transforms(bool inc_mirror) const
{
	vector<Transform>  ret;
	vector<int> hit_cache;
	
	vector<Vec3f> points = get_asym_unit_points(inc_mirror);
	// Warning, this is a gross hack because it is assuming that the asym_unit_points
	// returned by DSym are in a particular orientation with respect to symmetric axes
	// if the internals of DSym change it could change what we should do here...
	// but for the time being it will do
	if (inc_mirror && is_d_sym() && (get_nsym()/2 % 2 == 0)) {
		Dict delim = get_delimiters(false);
		float angle = (float)(EMConsts::deg2rad*float(delim["az_max"]));
		float y = -cos(angle);
		float x = sin(angle);
		points.push_back(Vec3f(x,y,0));
	}
	else if ( is_d_sym() && (get_nsym()/2 % 2 == 1)) {
		Dict delim = get_delimiters(false);
		float angle = float(delim["az_max"])/2.0f;
// 		cout << "Odd dsym using " << angle << endl;
		angle *= (float)EMConsts::deg2rad;
		float y = -cos(angle);
		float x = sin(angle);
		points.push_back(Vec3f(x,y,0));
		
		if ( inc_mirror ) {
			angle = 3.0f*(float(delim["az_max"]))/2.0f;
			angle *= (float)EMConsts::deg2rad;
			float y = -cos(angle);
			float x = sin(angle);
			points.push_back(Vec3f(x,y,0));
		}
	}
	
	typedef vector<Vec3f>::const_iterator const_point_it;
	for(const_point_it point = points.begin(); point != points.end(); ++point ) {
		
		for(int i = 1; i < get_nsym(); ++i) {
			
			if ( find(hit_cache.begin(),hit_cache.end(),i) != hit_cache.end() ) continue;
			Transform t = get_sym(i);
			Vec3f result = (*point)*t;
			
			if (is_platonic_sym()) {
				for(const_point_it tmp = points.begin(); tmp != points.end(); ++tmp ) {
					Vec3f tt = result-(*tmp);
					if (tt.squared_length() < 0.01) {
						hit_cache.push_back(i);
						ret.push_back(t);
					}
					
				}
			}
			else {
				result -= *point;
				if (result.squared_length() < 0.05) {
					hit_cache.push_back(i);
					ret.push_back(t);
				}
			}
		}
		
	}
	
	return ret;
}


// C Symmetry stuff 
Dict CSym::get_delimiters(const bool inc_mirror) const {
	Dict returnDict;
	// Get the parameters of interest
	int nsym = params.set_default("nsym",0);
	if ( nsym <= 0 ) throw InvalidValueException(nsym,"Error, you must specify a positive non zero nsym");
			
	if ( inc_mirror ) returnDict["alt_max"] = 180.0;
	else  returnDict["alt_max"] = 90.0;
		
	returnDict["az_max"] = 360.0/(float)nsym;
			
	return returnDict;
}

bool CSym::is_in_asym_unit(const float& altitude, const float& azimuth, const bool inc_mirror = false) const
{
	Dict d = get_delimiters(inc_mirror);
	float alt_max = d["alt_max"];
	float az_max = d["az_max"];
	
	int nsym = params.set_default("nsym",0);
	if ( nsym != 1 && azimuth < 0) return false;
	if ( altitude <= alt_max && azimuth <= az_max ) return true;
	return false;
}

vector<vector<Vec3f> > CSym::get_asym_unit_triangles(bool inc_mirror ) const{
	vector<Vec3f> v = get_asym_unit_points(inc_mirror);
	int nsym = params.set_default("nsym",0);
	
	vector<vector<Vec3f> > ret;
	if (v.size() == 0) return ret; // nsym == 1 and inc_mirror == true, this is the entire sphere!
	if (nsym == 1 && !inc_mirror) {
		Vec3f z(0,0,1);
		vector<Vec3f> tmp;
		tmp.push_back(z);
		tmp.push_back(v[1]);
		tmp.push_back(v[0]);
		ret.push_back(tmp);
		
		vector<Vec3f> tmp2;
		tmp2.push_back(z);
		tmp2.push_back(v[2]);
		tmp2.push_back(v[1]);
		ret.push_back(tmp2);
		
		vector<Vec3f> tmp3;
		tmp3.push_back(z);
		tmp3.push_back(v[3]);
		tmp3.push_back(v[2]);
		ret.push_back(tmp3);
		
		vector<Vec3f> tmp4;
		tmp4.push_back(z);
		tmp4.push_back(v[0]);
		tmp4.push_back(v[3]);
		ret.push_back(tmp4);
	}
	else if (nsym == 2 && inc_mirror) {
		Vec3f x(1,0,0);
		vector<Vec3f> tmp;
		tmp.push_back(v[1]);
		tmp.push_back(v[0]);
		tmp.push_back(x);
		ret.push_back(tmp);
		
		vector<Vec3f> tmp2;
		tmp2.push_back(v[2]);
		tmp2.push_back(v[1]);
		tmp2.push_back(x);
		ret.push_back(tmp2);
		
		vector<Vec3f> tmp3;
		tmp3.push_back(v[3]);
		tmp3.push_back(v[2]);
		tmp3.push_back(x);
		ret.push_back(tmp3);
		
		vector<Vec3f> tmp4;
		tmp4.push_back(v[0]);
		tmp4.push_back(v[3]);
		tmp4.push_back(x);
		ret.push_back(tmp4);
	}
	else if (nsym == 2 && !inc_mirror) {
		vector<Vec3f> tmp;
		tmp.push_back(v[0]);
		tmp.push_back(v[2]);
		tmp.push_back(v[1]);
		ret.push_back(tmp);
		
		vector<Vec3f> tmp2;
		tmp2.push_back(v[2]);
		tmp2.push_back(v[0]);
		tmp2.push_back(v[3]);
		ret.push_back(tmp2);
	}
	else if (v.size() == 3) {
			vector<Vec3f> tmp;
			tmp.push_back(v[0]);
			tmp.push_back(v[2]);
			tmp.push_back(v[1]);
			ret.push_back(tmp);
	}
	else if (v.size() == 4) {
		vector<Vec3f> tmp;
		tmp.push_back(v[0]);
		tmp.push_back(v[3]);
		tmp.push_back(v[1]);
		ret.push_back(tmp);
		
		vector<Vec3f> tmp2;
		tmp2.push_back(v[1]);
		tmp2.push_back(v[3]);
		tmp2.push_back(v[2]);
		ret.push_back(tmp2);
	}
	
	return ret;
}

vector<Vec3f> CSym::get_asym_unit_points(bool inc_mirror) const
{
	Dict delim = get_delimiters(inc_mirror);
	int nsym = params.set_default("nsym",0);
	vector<Vec3f> ret;
	
	if ( nsym == 1 ) {
		if (inc_mirror == false ) {
			ret.push_back(Vec3f(0,-1,0));
			ret.push_back(Vec3f(1,0,0));
			ret.push_back(Vec3f(0,1,0));
			ret.push_back(Vec3f(-1,0,0));
		}
		// else return ret; // an empty vector! this is fine
	}
	else if (nsym == 2 && !inc_mirror) {
		ret.push_back(Vec3f(0,0,1));
		ret.push_back(Vec3f(0,-1,0));
		ret.push_back(Vec3f(1,0,0));
		ret.push_back(Vec3f(0,1,0));
	}
	else {
		ret.push_back(Vec3f(0,0,1));
		ret.push_back(Vec3f(0,-1,0));
		if (inc_mirror == true) {
			ret.push_back(Vec3f(0,0,-1));
		}
		float angle = (float)(EMConsts::deg2rad*float(delim["az_max"]));
		float y = -cos(angle);
		float x = sin(angle);
		ret.push_back(Vec3f(x,y,0));
	}
	
	return ret;
	
}
		
Transform CSym::get_sym(const int n) const {
	int nsym = params.set_default("nsym",0);
	if ( nsym <= 0 ) throw InvalidValueException(n,"Error, you must specify a positive non zero nsym");
			
	Dict d("type","eman");
	// courtesy of Phil Baldwin
	d["az"] = (n%nsym) * 360.0f / nsym;
	d["alt"] = 0.0f;
	d["phi"] = 0.0f;
	return Transform(d);
}

// D symmetry stuff
Dict DSym::get_delimiters(const bool inc_mirror) const {
	Dict returnDict;
			
	// Get the parameters of interest
	int nsym = params.set_default("nsym",0);
	if ( nsym <= 0 ) throw InvalidValueException(nsym,"Error, you must specify a positive non zero nsym");
			
	returnDict["alt_max"] = 90.0;
		
	if ( inc_mirror )  returnDict["az_max"] = 360.0/(float)nsym;
	else returnDict["az_max"] = 180.0/(float)nsym;
				
	return returnDict;
}

bool DSym::is_in_asym_unit(const float& altitude, const float& azimuth, const bool inc_mirror = false) const
{
	Dict d = get_delimiters(inc_mirror);
	float alt_max = d["alt_max"];
	float az_max = d["az_max"];
	
	int nsym = params.set_default("nsym",0);
	
	if ( nsym == 1 && inc_mirror ) {
		 if (altitude >= 0 && altitude <= alt_max && azimuth <= az_max ) return true;
	}
	else {
		if ( altitude >= 0 && altitude <= alt_max && azimuth <= az_max && azimuth >= 0 ) return true;
	}
	return false;
}

Transform DSym::get_sym(const int n) const
{
	int nsym = 2*params.set_default("nsym",0);
	if ( nsym <= 0 ) throw InvalidValueException(n,"Error, you must specify a positive non zero nsym");
	
	Dict d("type","eman");
	// courtesy of Phil Baldwin
	if (n >= nsym / 2) {
		d["az"] = ( (n%nsym) - nsym/2) * 360.0f / (nsym / 2);
		d["alt"] = 180.0f;
		d["phi"] = 0.0f;
	}
	else {
		d["az"] = (n%nsym) * 360.0f / (nsym / 2);
		d["alt"] = 0.0f;
		d["phi"] = 0.0f;
	}
	return Transform(d);
}

vector<vector<Vec3f> > DSym::get_asym_unit_triangles(bool inc_mirror) const{
	vector<Vec3f> v = get_asym_unit_points(inc_mirror);
	int nsym = params.set_default("nsym",0);
	vector<vector<Vec3f> > ret;
	if ( (nsym == 1 && inc_mirror == false) || (nsym == 2 && inc_mirror)) {
			vector<Vec3f> tmp;
			tmp.push_back(v[0]);
			tmp.push_back(v[2]);
			tmp.push_back(v[1]);
			ret.push_back(tmp);
		
			vector<Vec3f> tmp2;
			tmp2.push_back(v[2]);
			tmp2.push_back(v[0]);
			tmp2.push_back(v[3]);
			ret.push_back(tmp2);
	}
	else if (nsym == 1) {
		Vec3f z(0,0,1);
		vector<Vec3f> tmp;
		tmp.push_back(z);
		tmp.push_back(v[1]);
		tmp.push_back(v[0]);
		ret.push_back(tmp);
	
		vector<Vec3f> tmp2;
		tmp2.push_back(z);
		tmp2.push_back(v[2]);
		tmp2.push_back(v[1]);
		ret.push_back(tmp2);
	
		vector<Vec3f> tmp3;
		tmp3.push_back(z);
		tmp3.push_back(v[3]);
		tmp3.push_back(v[2]);
		ret.push_back(tmp3);
	
		vector<Vec3f> tmp4;
		tmp4.push_back(z);
		tmp4.push_back(v[0]);
		tmp4.push_back(v[3]);
		ret.push_back(tmp4);
	}
	else {
// 		if v.size() == 3
		vector<Vec3f> tmp;
		tmp.push_back(v[0]);
		tmp.push_back(v[2]);
		tmp.push_back(v[1]);
		ret.push_back(tmp);
	}
	
	return ret;
}

vector<Vec3f> DSym::get_asym_unit_points(bool inc_mirror) const
{
	Dict delim = get_delimiters(inc_mirror);
	
	vector<Vec3f> ret;
	int nsym = params.set_default("nsym",0);
	if ( nsym == 1 ) {
		if (inc_mirror == false ) {
			ret.push_back(Vec3f(0,0,1));
			ret.push_back(Vec3f(0,-1,0));
			ret.push_back(Vec3f(1,0,0));
			ret.push_back(Vec3f(0,1,0));
		}
		else {
			ret.push_back(Vec3f(0,-1,0));
			ret.push_back(Vec3f(1,0,0));
			ret.push_back(Vec3f(0,1,0));
			ret.push_back(Vec3f(-1,0,0));
		}
	}
	else if ( nsym == 2 && inc_mirror ) {
		ret.push_back(Vec3f(0,0,1));
		ret.push_back(Vec3f(0,-1,0));
		ret.push_back(Vec3f(1,0,0));
		ret.push_back(Vec3f(0,1,0));
	}
	else {
		float angle = (float)(EMConsts::deg2rad*float(delim["az_max"]));
		ret.push_back(Vec3f(0,0,1));
		ret.push_back(Vec3f(0,-1,0));
		float y = -cos(angle);
		float x = sin(angle);
		ret.push_back(Vec3f(x,y,0));
	}
	
	return ret;
	
}

// H symmetry stuff
Dict HSym::get_delimiters(const bool) const {
	Dict returnDict;
	
	// Get the parameters of interest
	int nsym = params.set_default("nsym",0);
	if ( nsym <= 0 ) throw InvalidValueException(nsym,"Error, you must specify a positive non zero nsym");
	
	float equator_range = params.set_default("equator_range",5.0f);
	
	returnDict["alt_max"] = 90.0 + equator_range;
	
	returnDict["alt_min"] = 90.0;
		
	returnDict["az_max"] = 360.0/(float)nsym;

	return returnDict;
}

bool HSym::is_in_asym_unit(const float& altitude, const float& azimuth, const bool inc_mirror = false) const
{
	Dict d = get_delimiters(inc_mirror);
	float alt_max = d["alt_max"];
	float alt_min = d["alt_min"];
	
	if (inc_mirror) {
		float e = params.set_default("equator_range",5.0f);
		alt_min -= e;
	}
	
	float az_max = d["az_max"];
	
	if ( altitude >=alt_min && altitude <= alt_max && azimuth <= az_max && azimuth >= 0 ) return true;
	return false;
}

vector<vector<Vec3f> > HSym::get_asym_unit_triangles(bool inc_mirror) const{
	
	vector<vector<Vec3f> > ret;
	return ret;
}

vector<Vec3f> HSym::get_asym_unit_points(bool inc_mirror) const
{
	vector<Vec3f> ret;
	
	Dict delim = get_delimiters(inc_mirror);
	int nsym = params.set_default("nsym",0);
	float az = -(float)delim["az_max"];
	
	
	bool tracing_arcs = false;
	
	
	if ( !tracing_arcs) {
		Vec3f a(0,-1,0);
		ret.push_back(a);
		
		if ( nsym > 2 )	{
			Dict d("type","eman");
			d["phi"] = 0.0f;
			d["alt"] = 0.0f;
			d["az"] = az;
			Vec3f b = Transform(d)*a; 
			ret.push_back(b);
		}
		else
		{
			ret.push_back(Vec3f(1,0,0));
			
			ret.push_back(Vec3f(0,1,0));
			
			if ( nsym == 1 ) {
				ret.push_back(Vec3f(-1,0,0));
				ret.push_back(a);
			}
		}
	}
	return ret;
	
}

Transform HSym::get_sym(const int n) const
{
	float daz = params.set_default("daz",0.0f);
	float apix = params.set_default("apix",1.0f);
	float dz = params.set_default("tz", 0)/apix;
	
	Dict d("type","eman");
	// courtesy of Phil Baldwin
	d["az"] = n * daz;
	d["alt"] = 0.0f;
	d["phi"] = 0.0f;
	Transform ret(d);
	ret.set_trans(0,0,n*dz);
	return ret;
}

// Generic platonic symmetry stuff
void PlatonicSym::init()
{
	//See the manuscript "The Transform Class in Sparx and EMAN2", Baldwin & Penczek 2007. J. Struct. Biol. 157 (250-261)
	//In particular see pages 257-259
	//cap_sig is capital sigma in the Baldwin paper
	float cap_sig =  2.0f*M_PI/ get_max_csym();
	//In EMAN2 projection cap_sig is really az_max
	platonic_params["az_max"] = cap_sig;
			
	// Alpha is the angle between (immediately) neighborhing 3 fold axes of symmetry
	// This follows the conventions in the Baldwin paper
	float alpha = acos(1.0f/(sqrtf(3.0f)*tan(cap_sig/2.0f)));
	// In EMAN2 projection alpha is really al_maz
	platonic_params["alt_max"] = alpha;
			
	// This is half of "theta_c" as in the conventions of the Balwin paper. See also http://blake.bcm.edu/emanwiki/EMAN2/Symmetry.
	platonic_params["theta_c_on_two"] = 1.0/2.0*acos( cos(cap_sig)/(1.0-cos(cap_sig)));
	
}


Dict PlatonicSym::get_delimiters(const bool inc_mirror) const
{	
	Dict ret;
	ret["az_max"] = EMConsts::rad2deg * (float) platonic_params["az_max"];
	// For icos and oct symmetries, excluding the mirror means halving az_maz
	if ( inc_mirror == false ) 
		if ( get_name() ==  IcosahedralSym::NAME || get_name() == OctahedralSym::NAME ) 
			ret["az_max"] = 0.5*EMConsts::rad2deg * (float) platonic_params["az_max"];
		//else
		//the alt_max variable should probably be altered if the symmetry is tet, but
		//this is taken care of in TetSym::is_in_asym_unit
	
	ret["alt_max"] = EMConsts::rad2deg * (float) platonic_params["alt_max"];
	return ret;
}

//.Warning, this function only returns valid answers for octahedral and icosahedral symmetries.
bool PlatonicSym::is_in_asym_unit(const float& altitude, const float& azimuth, const bool inc_mirror) const
{
	Dict d = get_delimiters(inc_mirror);
	float alt_max = d["alt_max"];
	float az_max = d["az_max"];
	
	if ( altitude >= 0 &&  altitude <= alt_max && azimuth <= az_max && azimuth >= 0) {
	
		// Convert azimuth to radians
		float tmpaz = (float)(EMConsts::deg2rad * azimuth);
		
		float cap_sig = platonic_params["az_max"];
		float alt_max = platonic_params["alt_max"];
		if ( tmpaz > ( cap_sig/2.0 ) )tmpaz = cap_sig - tmpaz;
		
		float lower_alt_bound = platonic_alt_lower_bound(tmpaz, alt_max );
		
		// convert altitude to radians
		float tmpalt = (float)(EMConsts::deg2rad * altitude);
		if ( lower_alt_bound > tmpalt ) {
			if ( inc_mirror == false )
			{
				if ( cap_sig/2.0 < tmpaz ) return false;
				else return true;
			}
			else return true;
		}
		return false;
	}
	return false;
}

float PlatonicSym::platonic_alt_lower_bound(const float& azimuth, const float& alpha) const
{
	float cap_sig = platonic_params["az_max"];
	float theta_c_on_two = platonic_params["theta_c_on_two"];
	
	float baldwin_lower_alt_bound = sin(cap_sig/2.0f-azimuth)/tan(theta_c_on_two);
	baldwin_lower_alt_bound += sin(azimuth)/tan(alpha);
	baldwin_lower_alt_bound *= 1/sin(cap_sig/2.0f);
	baldwin_lower_alt_bound = atan(1/baldwin_lower_alt_bound);

	return baldwin_lower_alt_bound;
}

vector<vector<Vec3f> > PlatonicSym::get_asym_unit_triangles(bool inc_mirror) const{
	vector<Vec3f> v = get_asym_unit_points(inc_mirror);
	vector<vector<Vec3f> > ret;
	if (v.size() == 3) {
		vector<Vec3f> tmp;
		tmp.push_back(v[0]);
		tmp.push_back(v[2]);
		tmp.push_back(v[1]);
		ret.push_back(tmp);
	}
	else /* v.size() == 4*/ {
		vector<Vec3f> tmp;
		tmp.push_back(v[0]);
		tmp.push_back(v[2]);
		tmp.push_back(v[1]);
		ret.push_back(tmp);
		
		vector<Vec3f> tmp2;
		tmp2.push_back(v[0]);
		tmp2.push_back(v[3]);
		tmp2.push_back(v[2]);
		ret.push_back(tmp2);
	}
	
	return ret;
}

vector<Vec3f> PlatonicSym::get_asym_unit_points(bool inc_mirror) const
{
	vector<Vec3f> ret;
	
	Vec3f b = Vec3f(0,0,1);
	ret.push_back(b);
	float theta_c_on_two = (float)platonic_params["theta_c_on_two"]; // already in radians
	float theta_c = 2*theta_c_on_two;
	
	Vec3f c_on_two = Vec3f(0,-sin(theta_c_on_two),cos(theta_c_on_two));
	Vec3f c = Vec3f(0,-sin(theta_c),cos(theta_c));
	ret.push_back(c_on_two);
	
	float cap_sig = platonic_params["az_max"];
	Vec3f a = Vec3f(sin(theta_c)*sin(cap_sig),-sin(theta_c)*cos(cap_sig),cos(theta_c));
	
	Vec3f f = a+b+c;
	f.normalize();
		
	ret.push_back(f);
	
	if ( inc_mirror ) {
		Vec3f a_on_two = Vec3f(sin(theta_c_on_two)*sin(cap_sig),-sin(theta_c_on_two)*cos(cap_sig),cos(theta_c_on_two));
		ret.push_back(a_on_two);
	}
	
	if ( get_az_alignment_offset() != 0 ) {
		Dict d("type","eman");
		d["az"] = get_az_alignment_offset();
		d["phi"] = 0.0f;
		d["alt"] = 0.0f;
		Transform t(d);
		for (vector<Vec3f>::iterator it = ret.begin(); it != ret.end(); ++it )
		{
			*it = (*it)*t;
		} 
	}
// 	
	return ret;
	
}

float IcosahedralSym::get_az_alignment_offset() const { return 234.0; } // This offset positions a 3 fold axis on the positive x axis

Transform IcosahedralSym::get_sym(const int n) const
{
	// These rotations courtesy of Phil Baldwin
	static double  lvl0=0.; //  there is one pentagon on top; five-fold along z
	static double  lvl1= 63.4349; // that is atan(2)  // there are 5 pentagons with centers at this height (angle)
	static double  lvl2=116.5651; //that is 180-lvl1  // there are 5 pentagons with centers at this height (angle)
	static double lvl3=180.0;
	
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
	
	int idx = n % 60;
	Dict d("type","eman");
// 	Transform3D ret;
	if (get_az_alignment_offset() == 234.0) {
		d["az"] =(float)ICOS[idx * 3 ]+90;
		d["alt"] = (float)ICOS[idx * 3 + 1];
		d["phi"] = (float)ICOS[idx * 3 + 2]-90;
// 		ret.set_rotation((float)ICOS[idx * 3 ]+90,(float)ICOS[idx * 3 + 1], (float)ICOS[idx * 3 + 2]-90);
	}
	else {
		d["az"] =(float)(float)ICOS[idx * 3 ];
		d["alt"] = (float)ICOS[idx * 3 + 1];
		d["phi"] = (float)ICOS[idx * 3 + 2];
// 		ret.set_rotation((float)ICOS[idx * 3 ],(float)ICOS[idx * 3 + 1], (float)ICOS[idx * 3 + 2]);
	}
	
// 	ret.set_rotation((float)ICOS[idx * 3 ],(float)ICOS[idx * 3 + 1], (float)ICOS[idx * 3 + 2]);
// 	if ( get_az_alignment_offset() != 0 ) {
// 		Transform3D t(get_az_alignment_offset(),0,0);
// 		ret = t*ret;
// 	}
	return Transform(d);
	
}

Transform OctahedralSym::get_sym(const int n) const
{
	// These rotations courtesy of Phil Baldwin
	// We have placed the OCT symmetry group with a face along the z-axis
	static double lvl0=0.;
	static double lvl1=90.;
	static double lvl2=180.;

	static double OCT[72] = {// This is with a face of a cube along z 
		0,lvl0,0,   0,lvl0,90,    0,lvl0,180,    0,lvl0,270,
		0,lvl1,0,   0,lvl1,90,    0,lvl1,180,    0,lvl1,270,
		90,lvl1,0,  90,lvl1,90,   90,lvl1,180,   90,lvl1,270,
		180,lvl1,0, 180,lvl1,90,  180,lvl1,180,  180,lvl1,270,
		270,lvl1,0, 270,lvl1,90,  270,lvl1,180,  270,lvl1,270,
		0,lvl2,0,   0,lvl2,90,    0,lvl2,180,    0,lvl2,270
	};
	
	int idx = n % 24;
// 	Transform3D ret;
// 	ret.set_rotation((float)OCT[idx * 3 ],(float)OCT[idx * 3 + 1], (float)OCT[idx * 3 + 2] );
	Dict d("type","eman");
	d["az"] = (float)OCT[idx * 3 ];
	d["alt"] = (float)OCT[idx * 3 + 1];
	d["phi"] = (float)OCT[idx * 3 + 2];
	return Transform(d);
	
}

float TetrahedralSym::get_az_alignment_offset() const { return  0.0; }

bool TetrahedralSym::is_in_asym_unit(const float& altitude, const float& azimuth, const bool inc_mirror) const
{
	Dict d = get_delimiters(inc_mirror);
	float alt_max = d["alt_max"];
	float az_max = d["az_max"];
	
	if ( altitude >= 0 &&  altitude <= alt_max && azimuth <= az_max && azimuth >= 0) {
		// convert azimuth to radians
		float tmpaz = (float)(EMConsts::deg2rad * azimuth);
			
		float cap_sig = platonic_params["az_max"];
		float alt_max = platonic_params["alt_max"];
		if ( tmpaz > ( cap_sig/2.0 ) )tmpaz = cap_sig - tmpaz;
		
		float lower_alt_bound = platonic_alt_lower_bound(tmpaz, alt_max );
		
		// convert altitude to radians
		float tmpalt = (float)(EMConsts::deg2rad * altitude);
		if ( lower_alt_bound > tmpalt ) {
			if ( !inc_mirror ) {
				float upper_alt_bound = platonic_alt_lower_bound( tmpaz, alt_max/2.0f);
				// you could change the "<" to a ">" here to get the other mirror part of the asym unit
				if ( upper_alt_bound < tmpalt ) return false;
				else return true;
			}
			else return true;
		}
		return false;
	}
	else return false;
}


Transform TetrahedralSym::get_sym(const int n) const
{
	// These rotations courtesy of Phil Baldwin
	 // It has n=m=3; F=4, E=6=nF/2, V=4=nF/m
	static double lvl0=0;         // There is a face along z
	static double lvl1=109.4712;  //  that is acos(-1/3)  // There  are 3 faces at this angle

	static double TET[36] = {// This is with the face along z 
	      0,lvl0,0,   0,lvl0,120,    0,lvl0,240,
	      0,lvl1,60,   0,lvl1,180,    0,lvl1,300,
	    120,lvl1,60, 120,lvl1,180,  120,lvl1,300,
	    240,lvl1,60, 240,lvl1,180,  240,lvl1,300
	};
// 	
	int idx = n % 12;
// 	Transform3D ret;
// 	ret.set_rotation((float)TET[idx * 3 ],(float)TET[idx * 3 + 1], (float)TET[idx * 3 + 2] );
// 	return ret;
	
	Dict d("type","eman");
	d["az"] = (float)TET[idx * 3 ];
	d["alt"] = (float)TET[idx * 3 + 1];
	d["phi"] = (float)TET[idx * 3 + 2];
	return Transform(d);
	
}


vector<Vec3f> TetrahedralSym::get_asym_unit_points(bool inc_mirror) const
{
	vector<Vec3f> ret;
	
	Vec3f b = Vec3f(0,0,1);
	ret.push_back(b);
	float theta_c_on_two = (float)platonic_params["theta_c_on_two"]; // already in radians
	float theta_c = 2*theta_c_on_two;
	
	Vec3f c_on_two = Vec3f(0,-sin(theta_c_on_two),cos(theta_c_on_two));
	Vec3f c = Vec3f(0,-sin(theta_c),cos(theta_c));
	ret.push_back(c_on_two);
	float cap_sig = platonic_params["az_max"];
	if ( inc_mirror ) {
		Vec3f a = Vec3f(sin(theta_c)*sin(cap_sig),-sin(theta_c)*cos(cap_sig),cos(theta_c));
		
		Vec3f f = a+b+c;
		f.normalize();
			
		ret.push_back(f);
	}
	
	Vec3f a_on_two = Vec3f(sin(theta_c_on_two)*sin(cap_sig),-sin(theta_c_on_two)*cos(cap_sig),cos(theta_c_on_two));
	ret.push_back(a_on_two);

	
	if ( get_az_alignment_offset() != 0 ) {
		Dict d("type","eman");
		d["az"] = get_az_alignment_offset();
		d["phi"] = 0.0f;
		d["alt"] = 0.0f;
		Transform t(d);
		for (vector<Vec3f>::iterator it = ret.begin(); it != ret.end(); ++it )
		{
			*it = (*it)*t;
		} 
	}
	
	return ret;
	
}

/* vim: set ts=4 noet: */
