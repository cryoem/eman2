/*
 * Author: David Woolford, 09/23/2008 (woolford@bcm.edu)
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
 * source code. Additional authorship citations may be added, but existingx
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

#include "symmetry.h"
#include "transform.h"
#include "vec3.h"
#include "exception.h"
#include "util.h"

using namespace EMAN;

const string CSym::NAME = "c";
const string DSym::NAME = "d";
const string HSym::NAME = "h";
const string TetrahedralSym::NAME = "tet";
const string OctahedralSym::NAME = "oct";
const string IcosahedralSym::NAME = "icos";
const string Icosahedral2Sym::NAME = "icos2";
const string EmanOrientationGenerator::NAME = "eman";
const string SaffOrientationGenerator::NAME = "saff";
const string EvenOrientationGenerator::NAME = "even";
const string RandomOrientationGenerator::NAME = "rand";
const string OptimumOrientationGenerator::NAME = "opt";

template <> Factory < Symmetry3D >::Factory()
{
	force_add<CSym>();
	force_add<DSym>();
	force_add<HSym>();
	force_add<TetrahedralSym>();
	force_add<OctahedralSym>();
	force_add<IcosahedralSym>();
	force_add<Icosahedral2Sym>();
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
Symmetry3D* Factory < Symmetry3D >::get(const string & instancename_)
{
	init();

	string instancename = Util::str_to_lower(instancename_);

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
			int nstart=1,nsym=1;
			float daz=5.,tz=5.,maxtilt=1.0;
			string temp;
			temp=instancename;
			temp.erase(0,1);
			if (sscanf(temp.c_str(),"%d:%d:%f:%f:%f",&nsym,&nstart,&daz,&tz,&maxtilt)<4) {
				sscanf(temp.c_str(),"%d,%d,%f,%f,%f",&nsym,&nstart,&daz,&tz,&maxtilt);
			}
			parms["nstart"]=nstart;
			parms["nsym"]=nsym;
			parms["daz"]=daz;
			parms["tz"]=tz;
			parms["maxtilt"]=maxtilt;
			return get("h",parms);
		}

// 		delete lc;
	}
	else if ( instancename == "icos" || instancename == "oct" || instancename == "tet" || instancename == "icos2" || instancename == "icos5" )
	{
		if (instancename == "icos5") {
			instancename = "icos";
		}

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
	force_add<EmanOrientationGenerator>();
	force_add<RandomOrientationGenerator>();
	force_add<EvenOrientationGenerator>();
	force_add<SaffOrientationGenerator>();
	force_add<OptimumOrientationGenerator>();
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
	OrientationGenerator *g = Factory < OrientationGenerator >::get(Util::str_to_lower(generatorname), parms);
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

void OrientationGenerator::get_az_max(const Symmetry3D* const sym, const float& altmax, const bool inc_mirror, const float& alt_iterator,const float& h,bool& d_odd_mirror_flag, float& azmax_adjusted) const
{

	if ( sym->is_d_sym() && alt_iterator == altmax && ( (sym->get_nsym())/2 % 2 == 1 )) {
		if (inc_mirror) {
			azmax_adjusted /= 4.0f;
			d_odd_mirror_flag = true;
		}
		else azmax_adjusted /= 2.0f;
	}
	else if (sym->is_d_sym() && alt_iterator == altmax && ( (sym->get_nsym())/2 % 2 == 0 ) && inc_mirror) {
		azmax_adjusted /= 2.0f;
	}
	// if this is odd c symmetry, and we're at the equator, and we're excluding the mirror then
	// half the equator is redundant (it is the mirror of the other half)
	else if (sym->is_c_sym() && !inc_mirror && alt_iterator == altmax && (sym->get_nsym() % 2 == 1 ) ){
		azmax_adjusted /= 2.0f;
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

}


//vector<Vec3f> OrientationGenerator::get_graph_opt(const Symmetry3D* const sym) const {
//	bool inc_mirror = params.set_default("inc_mirror",false);
//	vector<Vec3f> seeds = get_asym_unit_points(inc_mirror);
//	vector<Vec3i> connections;
//	if (seeds.size() == 3) {
//		connections.push_back(Vec2i(0,1,2));
//	} else {
//		throw;
//
//	}
//}


float OrientationGenerator::get_optimal_delta(const Symmetry3D* const sym, const int& n) const
{

//	float delta_soln = 360.0f/sym->get_max_csym();
	float delta_soln = 180.0f;
	float delta_upper_bound = delta_soln;
	float delta_lower_bound = 0.0f;

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
	float phi = 0.0f;
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
			for ( float p = phitoo; p <= 360.0f-phitoo; p+= phitoo )
			{
				d["phi"] = fmod(phi+p,360);
				Transform t(d);
				v.push_back(t);
			}
		}
	}
	return true;
}

float EmanOrientationGenerator::get_az_delta(const float& delta,const float& altitude, const int) const
{
	// convert altitude into radians
	float tmp = (float)(EMConsts::deg2rad * altitude);

	// This is taken from EMAN1 project3d.C
	// This wasn't working like it was supposed to. Rather than
	// figuring it out, I'm just replacing it --steve
/*	float h=floor(360.0f/(delta*1.1547f));	// the 1.1547 makes the overall distribution more like a hexagonal mesh
	h=(int)floor(h*sin(tmp)+.5f);
	if (h==0) h=1;
	h=abs(maxcsym)*floor(h/(float)abs(maxcsym)+.5f);
	if ( h == 0 ) h = (float)maxcsym;
	h=2.0f*M_PI/h;

	return (float)(EMConsts::rad2deg*h);*/

	return altitude==0?360.0f:delta/sin(tmp);

}


int EmanOrientationGenerator::get_orientations_tally(const Symmetry3D* const sym, const float& delta) const
{
	//FIXME THIS IS SO SIMILAR TO THE gen_orientations function that they should be probably use
	// a common routine - SAME ISSUE FOR OTHER ORIENTATION GENERATORS
	bool inc_mirror = params.set_default("inc_mirror",false);
	bool breaksym = params.set_default("breaksym",false);
	Dict delimiters = sym->get_delimiters(inc_mirror);
	float altmax = delimiters["alt_max"];
	float azmax = delimiters["az_max"];

	float paltmin = params.set_default("alt_min",0.0f);
	float paltmax = params.set_default("alt_max",180.0f);
	if (altmax>paltmax) altmax=paltmax;
	
	float alt_iterator = 0.0f;

	// #If it's a h symmetry then the alt iterator starts at very close
	// #to the altmax... the object is a h symmetry then it knows its alt_min...
	if (sym->is_h_sym()) alt_iterator = delimiters["alt_min"];

	int tally = 0;
	while ( alt_iterator <= altmax ) {
		float h = get_az_delta(delta,alt_iterator, sym->get_max_csym() );

		// not sure what this does code taken from EMAN1 - FIXME original author add comments
		if ( (alt_iterator > 0) && ( (azmax/h) < 2.8f) ) h = azmax / 2.1f;
		else if (alt_iterator == 0) h = azmax;

		float az_iterator = 0.0f;

		float azmax_adjusted = azmax;
		bool d_odd_mirror_flag = false;
		get_az_max(sym,altmax, inc_mirror,alt_iterator, h,d_odd_mirror_flag, azmax_adjusted);
		if (alt_iterator<paltmin) { alt_iterator += delta; continue; }
		
		while ( az_iterator <= azmax_adjusted ) {
			// FIXME: add an intelligent comment - this was copied from old code
//			if ( az_iterator > 180.0 && alt_iterator > 180.0/(2.0-0.001) && alt_iterator < 180.0/(2.0+0.001) ) {
//				az_iterator +=  h;
//				continue;
//			}

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
			if ( (az_iterator > azmax_adjusted) && d_odd_mirror_flag) {
				azmax_adjusted = azmax;
				az_iterator += azmax/2.0f;
			}
		}
		alt_iterator += delta;
	}
	
	if (breaksym) return tally*sym->get_nsym();
	return tally;
}

vector<Transform> EmanOrientationGenerator::gen_orientations(const Symmetry3D* const sym) const
{
	float delta = params.set_default("delta", 0.0f);
	int n = params.set_default("n", 0);
	bool breaksym = params.set_default("breaksym",false);

	if ( delta <= 0 && n <= 0 ) throw InvalidParameterException("Error, you must specify a positive non-zero delta or n");
	if ( delta > 0 && n > 0 ) throw InvalidParameterException("Error, the delta and the n arguments are mutually exclusive");

	if ( n > 0 ) {
		delta = get_optimal_delta(sym,n);
	}

	bool inc_mirror = params.set_default("inc_mirror",false);
	bool inc_mirror_real = inc_mirror;
	if (breaksym) inc_mirror=true;		// we need to enable mirror generation, then strip them out at the end, or things don't work right...
	Dict delimiters = sym->get_delimiters(inc_mirror);
	float altmax = delimiters["alt_max"];
	float azmax = delimiters["az_max"];

	float paltmin = params.set_default("alt_min",0.0f);
	float paltmax = params.set_default("alt_max",180.0f);
	if (altmax>paltmax) altmax=paltmax;

	bool perturb = params.set_default("perturb",true);

	float alt_iterator = 0.0f;

	// #If it's a h symmetry then the alt iterator starts at very close
	// #to the altmax... the object is a h symmetry then it knows its alt_min...
	if (sym->is_h_sym()) alt_iterator = delimiters["alt_min"];

	vector<Transform> ret;
	while ( alt_iterator <= altmax ) {
		float h = get_az_delta(delta,alt_iterator, sym->get_max_csym() );

		// not sure what this does code taken from EMAN1 - FIXME original author add comments
		if ( (alt_iterator > 0) && ( (azmax/h) < 2.8f) ) h = azmax / 2.1f;
		else if (alt_iterator == 0) h = azmax;

		float az_iterator = 0.0f;

		float azmax_adjusted = azmax;

		bool d_odd_mirror_flag = false;
		get_az_max(sym,altmax, inc_mirror,alt_iterator, h,d_odd_mirror_flag, azmax_adjusted);
		if (alt_iterator<paltmin) { alt_iterator += delta; continue; }


		while ( az_iterator <= azmax_adjusted ) {
			// FIXME: add an intelligent comment - this was copied from old code
// 			if ( az_iterator > 180.0 && alt_iterator > 180.0/(2.0-0.001) && alt_iterator < 180.0/(2.0+0.001) ) {
// 				az_iterator +=  h;
// 				continue;
// 			}
// 			// Now that I am handling the boundaries very specifically, I don't think we need
			// the above if statement. But I am leaving it there in case I need to reconsider.

			if (alt_iterator == 0 && az_iterator > 0){
				az_iterator += h;
				continue; // We only ever need to generate on orientation at alt=0
			}


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
//printf("%f %f/n",alt_soln,az_soln);
			if ( perturb &&  alt_soln != 0 ) {
				alt_soln += Util::get_gauss_rand(0.0f,.125f*delta);
				az_soln += Util::get_gauss_rand(0.0f,h*.125f);
			}

			add_orientation(ret,az_soln,alt_soln);

			// Add helical symmetry orientations on the other side of the equator (if we're including
			// mirror orientations)
			if ( sym->is_h_sym() && inc_mirror && alt_iterator != (float) delimiters["alt_min"] ) {
				add_orientation(ret, az_soln,2.0f*(float)delimiters["alt_min"]-alt_soln);
			}

			az_iterator += h;
			if ( (az_iterator > azmax_adjusted) && d_odd_mirror_flag) {
				azmax_adjusted = azmax;
				az_iterator += azmax/2.0f;
			}
		}
		alt_iterator += delta;
	}
	
	// With breaksym, values are generated for one asymmetric unit as if symmetry were imposed, then
	// the symmetry equivalent points are generated. Used with asymmetric structures with pseudosymmetry
	if (breaksym) {
		// no iterators here since we are making the list longer as we go
		int nwithsym=ret.size();	// transforms in one asym unit
		int nsym=sym->get_nsym();	// number of asymmetric units to generate
		for (int j=1; j<nsym; j++) {
			Transform t=sym->get_sym(j);
			for (int i=0; i<nwithsym; i++) {
				ret.push_back(ret[i]*t);		// add the symmetry modified transform to the end of the vector
			}
		}
		
		// Now we get rid of anything in the bottom half of the unit sphere if requested
		if (!inc_mirror_real) {
			vector<Transform> ret2;
			for (vector<Transform>::iterator t=ret.begin(); t!=ret.end(); ++t) {
				if ((*t)[2][2]>=0) ret2.push_back(*t); 
//				printf("%f\n",t[2][2]);
			}
			return ret2;
		}
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
		if ( s > 1.0f ) continue;
		float alpha = 2.0f*sqrtf(1.0f-s);
		float x = alpha * u1;
		float y = alpha * u2;
		float z = 2.0f*s-1.0f;

		float altitude = (float)EMConsts::rad2deg*acos(z);
		float azimuth = (float)EMConsts::rad2deg*atan2(y,x);

		float phi = 0.0f;
		if ( phitoo ) phi = Util::get_frand(0.0f,359.9999f);

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

	float altmin = 0.0f;
	// #If it's a h symmetry then the alt iterator starts at very close
	// #to the altmax... the object is a h symmetry then it knows its alt_min...
	if (sym->is_h_sym()) altmin = delimiters["alt_min"];

	int tally = 0;

	for (float alt = altmin; alt <= altmax; alt += delta) {
		float detaz;
		int lt;
		if ((0.0f == alt)||(180.0f == alt)) {
			detaz = 360.0f;
			lt = 1;
		} else {
			detaz = delta/(float)sin(alt*EMConsts::deg2rad);
			lt = int(azmax/detaz)-1;
			if (lt < 1) lt = 1;
			detaz = azmax/(float)lt;
		}
//		bool d_odd_mirror_flag = false;
//		get_az_max(sym,altmax, inc_mirror,alt, lt,d_odd_mirror_flag, detaz);
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

	float altmin = 0.0f;
	// If it's a h symmetry then the alt iterator starts at very close
	// to the altmax... the object is a h symmetry then it knows its alt_min...
	if (sym->is_h_sym()) altmin = delimiters["alt_min"];

	vector<Transform> ret;

	for (float alt = altmin; alt <= altmax; alt += delta) {
		float detaz;
		int lt;
		if ((0.0f == alt)||(180.0f == alt)) {
			detaz = 360.0f;
			lt = 1;
		} else {
			detaz = delta/(float)sin(alt*EMConsts::deg2rad);
			lt = int(azmax/detaz)-1;
			if (lt < 1) lt = 1;
			detaz = azmax/(float)lt;
		}
//		bool d_odd_mirror_flag = false;
//		get_az_max(sym,altmax, inc_mirror,alt, lt,d_odd_mirror_flag, detaz);

		for (int i = 0; i < lt; i++) {
			float az = (float)i*detaz;
			if (sym->is_platonic_sym()) {
				if ( sym->is_in_asym_unit(alt, az,inc_mirror) == false ) continue;
				else {
					az += sym->get_az_alignment_offset(); // Align to the symmetry axes
				}
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

	float altmin = 0.0f;
	// #If it's a h symmetry then the alt iterator starts at very close
	// #to the altmax... the object is a h symmetry then it knows its alt_min...
	if (sym->is_h_sym()){
		altmin = delimiters["alt_min"];
		if (inc_mirror) {
			altmin -= (float) sym->get_params()["maxtilt"];
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

	float altmin = 0.0f;
	// #If it's a h symmetry then the alt iterator starts at very close
	// #to the altmax... the object is a h symmetry then it knows its alt_min...
	if (sym->is_h_sym()){
		altmin = delimiters["alt_min"];
		if (inc_mirror) {
			altmin -= (float) sym->get_params()["maxtilt"];
		}
	}

	float Deltaz = (float)(cos(altmax*EMConsts::deg2rad)-cos(altmin*EMConsts::deg2rad));
	float s = delta*M_PI/180.0f;
	float NFactor = 3.6f/s;
	float wedgeFactor = fabs( Deltaz*(azmax)/720.0f) ;
	int NumPoints   =  static_cast<int> (NFactor*NFactor*wedgeFactor);

	vector<Transform> ret;

	if (!sym->is_h_sym()) add_orientation(ret,0,0);
	float az = 0.0f;
	float dz = (float)cos(altmin*EMConsts::deg2rad);
	for(int i = 1; i < NumPoints; ++i ){
		float z = dz + Deltaz* (float)i/ float(NumPoints-1);
		float r= sqrt(1.0f-z*z);
		az = fmod(az + delta/r,azmax);
		float alt = (float)(acos(z)*EMConsts::rad2deg);
		if (sym->is_platonic_sym()) {
			if ( sym->is_in_asym_unit(alt,az,inc_mirror) == false ) continue;
			else {
				az += sym->get_az_alignment_offset(); // Align to the symmetry axes
			}
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
		int max_it = 100;
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



Symmetry3D::Symmetry3D() : cached_au_planes(0),cache_size(0),num_triangles(0),au_sym_triangles() {}
Symmetry3D::~Symmetry3D() {
	if (cached_au_planes != 0 ) {
		delete_au_planes();
	}
}

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
//		throw;
		return t;
	}

	// Get the symmetry operation corresponding to the intersection asymmetric unit
	Transform nt = get_sym(soln);
	// Transpose it (invert it)
	nt.invert();
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
	o.invert();
	// Figure out where the point would end up. No we could just as easily not transpose and do
	// left multiplation (as in what occurs in the FourierReconstructor during slice insertion)
	p = o*p;

	return point_in_which_asym_unit(p);
}


void Symmetry3D::cache_au_planes() const {
	if (cached_au_planes == 0 ) {
		vector< vector<Vec3f> > au_triangles = get_asym_unit_triangles(true);
		num_triangles = au_triangles.size();
		cache_size = get_nsym()*au_triangles.size();

		cached_au_planes = new float*[cache_size];
		float** fit = cached_au_planes;
		for(int i =0; i < cache_size; ++i,++fit) {
			float *t = new float[4];
			*fit = t;
		}


		int k = 0;
		for(int i = 0; i < get_nsym(); ++i) {

			for( ncit it = au_triangles.begin(); it != au_triangles.end(); ++it, ++k)
			{
				// For each given triangle
				vector<Vec3f> points = *it;
				if ( i != 0 ) {
					for (vector<Vec3f>::iterator iit = points.begin(); iit != points.end(); ++iit ) {
						// Rotate the points in the triangle so that the triangle occupies the
						// space of the current asymmetric unit
						*iit = (*iit)*get_sym(i);
					}
				}

				au_sym_triangles.push_back(points);

				// Determine the equation of the plane for the points, store it in plane
				Util::equation_of_plane(points[0],points[2],points[1],cached_au_planes[k]);
			}
		}
	}
	else throw UnexpectedBehaviorException("Attempt to generate a cache when cache exists");
}

void Symmetry3D::delete_au_planes() {
	if (cached_au_planes == 0 ) throw UnexpectedBehaviorException("Attempt to delete a cache that does not exist");
	float** fit = cached_au_planes;
	for(int i =0; i < cache_size; ++i,++fit) {
		if (*fit == 0) throw UnexpectedBehaviorException("Attempt to delete a cache that does not exist");
		delete [] *fit;
		*fit = 0;
	}

	delete [] cached_au_planes;
	cached_au_planes = 0;
}

int Symmetry3D::point_in_which_asym_unit(const Vec3f& p) const
{
	if (cached_au_planes == 0) {
		cache_au_planes();
	}
	
	float epsNow=0.01f;
	int k = 0;
	for(int i = 0; i < get_nsym(); ++i) {
		for( int j = 0; j < num_triangles; ++j,++k) {
			vector<Vec3f> points = au_sym_triangles[k];

			float* plane = cached_au_planes[k];
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
			if ( s >= -epsNow && t >= -epsNow && (s+t) <= 1+epsNow ) {
//				cout << " i " << i << " j " << j << " s " << s  << " t " << t << " s+t " << s+t << endl;
				return i;
			}
		}
	}

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
					if (tt.squared_length() < 0.01f) {
						hit_cache.push_back(i);
						ret.push_back(t);
					}

				}
			}
			else {
				result -= *point;
				if (result.squared_length() < 0.05f) {
					hit_cache.push_back(i);
					ret.push_back(t);
				}
			}
		}

	}

	return ret;
}


vector<Transform> Symmetry3D::get_syms() const
{

	vector<Transform> ret;
// 	if (t.is_identity()) {
		for(int i = 0; i < get_nsym(); ++i ) {
			ret.push_back(get_sym(i));
		}
// 	} else {
// 		for(int i = 0; i < get_nsym(); ++i ) {
// 			ret.push_back(get_sym(i)*t);
// 		}
// 	}
	return ret;
}

vector<Transform> Symmetry3D::get_symmetries(const string& symmetry)
{
	Symmetry3D* sym = Factory<Symmetry3D>::get(Util::str_to_lower(symmetry));
	vector<Transform> ret = sym->get_syms();
	delete sym;
	return ret;
}

// C Symmetry stuff
Dict CSym::get_delimiters(const bool inc_mirror) const {
	Dict returnDict;
	// Get the parameters of interest
	int nsym = params.set_default("nsym",0);
	if ( nsym <= 0 ) throw InvalidValueException(nsym,"Error, you must specify a positive non zero nsym");

	if ( inc_mirror ) returnDict["alt_max"] = 180.0f;
	else  returnDict["alt_max"] = 90.0f;

	returnDict["az_max"] = 360.0f/(float)nsym;

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

	returnDict["alt_max"] = 90.0f;

	if ( inc_mirror )  returnDict["az_max"] = 360.0f/(float)nsym;
	else returnDict["az_max"] = 180.0f/(float)nsym;

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

	float maxtilt = params.set_default("maxtilt",5.0f);

	returnDict["alt_max"] = 90.0f;

	returnDict["alt_min"] = 90.0f - maxtilt;

	returnDict["az_max"] = 360.0f;

	return returnDict;
}

bool HSym::is_in_asym_unit(const float& altitude, const float& azimuth, const bool inc_mirror = false) const
{
	Dict d = get_delimiters(inc_mirror);
	float alt_max = d["alt_max"];
	float alt_min = d["alt_min"];

	if (inc_mirror) {
		float e = params.set_default("maxtilt",5.0f);
		alt_min -= e;
	}

	float az_max = d["az_max"];

	if ( altitude >=alt_min && altitude <= alt_max && azimuth <= az_max && azimuth >= 0 ) return true;
	return false;
}

vector<vector<Vec3f> > HSym::get_asym_unit_triangles(bool ) const{

	vector<vector<Vec3f> > ret;
	return ret;
}

vector<Vec3f> HSym::get_asym_unit_points(bool inc_mirror) const
{
	vector<Vec3f> ret;

	Dict delim = get_delimiters(inc_mirror);
	int nsym = params.set_default("nsym",1);
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
	int nstart=params["nstart"];
	//int nsym=params["nsym"];
	float apix = params.set_default("apix",1.0f);
	float daz= params["daz"];
	float tz=params["tz"];
	float dz=tz/apix;
	Dict d("type","eman");

	// courtesy of Phil Baldwin
	//d["az"] = (n%nsym) * 360.0f / nsym;
	//d["az"]=(((int) n/hsym)%nstart)*360.f/nstart+(n%hsym)*daz;
	//d["az"] = n * daz;
	d["az"]=(n%nstart)*(360.0/nstart)+floor(float(n)/nstart)*daz;	// corrected by steve, 7/21/11. No dependency on nsym
	d["alt"] = 0.0f;
	d["phi"] = 0.0f;
	Transform ret(d);
	ret.set_trans(0,0,(n/nstart)*dz);
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
	platonic_params["theta_c_on_two"] = 1.0f/2.0f*acos( cos(cap_sig)/(1.0f-cos(cap_sig)));

}


Dict PlatonicSym::get_delimiters(const bool inc_mirror) const
{
	Dict ret;
	ret["az_max"] = EMConsts::rad2deg * (float) platonic_params["az_max"];
	// For icos and oct symmetries, excluding the mirror means halving az_maz
	if ( inc_mirror == false )
		if ( get_name() ==  IcosahedralSym::NAME || get_name() == OctahedralSym::NAME )
			ret["az_max"] = 0.5f*EMConsts::rad2deg * (float) platonic_params["az_max"];
		//else
		//the alt_max variable should probably be altered if the symmetry is tet, but
		//this is taken care of in TetSym::is_in_asym_unit

	ret["alt_max"] = (float)(EMConsts::rad2deg * (float) platonic_params["alt_max"]);
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
		if ( tmpaz > ( cap_sig/2.0f ) )tmpaz = cap_sig - tmpaz;

		float lower_alt_bound = platonic_alt_lower_bound(tmpaz, alt_max );

		// convert altitude to radians
		float tmpalt = (float)(EMConsts::deg2rad * altitude);
		if ( lower_alt_bound > tmpalt ) {
			if ( inc_mirror == false )
			{
				if ( cap_sig/2.0f < tmpaz ) return false;
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

float Icosahedral2Sym::get_az_alignment_offset() const { return 234.0; } // This offset positions a 3 fold axis on the positive x axis (??? copied from IcosahedralSym)

Transform Icosahedral2Sym::get_sym(const int n) const
{
	static float matrices[60*9] = {
		1, 0, 0, 0, 1, 0, 0, 0, 1,
		0.30902, -0.80902, 0.5, 0.80902, 0.5, 0.30902, -0.5, 0.30902, 0.80902,
		-0.80902, -0.5, 0.30902, 0.5, -0.30902, 0.80902, -0.30902, 0.80902, 0.5,
		-0.80902, 0.5, -0.30902, -0.5, -0.30902, 0.80902, 0.30902, 0.80902, 0.5,
		0.30902, 0.80902, -0.5, -0.80902, 0.5, 0.30902, 0.5, 0.30902, 0.80902,
		-1, 0, 0, 0, -1, 0, 0, 0, 1,
		-0.30902, -0.80902, 0.5, 0.80902, -0.5, -0.30902, 0.5, 0.30902, 0.80902,
		0.30902, 0.80902, 0.5, -0.80902, 0.5, -0.30902, -0.5, -0.30902, 0.80902,
		-0.30902, 0.80902, 0.5, -0.80902, -0.5, 0.30902, 0.5, -0.30902, 0.80902,
		-0.5, 0.30902, 0.80902, 0.30902, -0.80902, 0.5, 0.80902, 0.5, 0.30902,
		0.5, -0.30902, 0.80902, -0.30902, 0.80902, 0.5, -0.80902, -0.5, 0.30902,
		0.80902, 0.5, 0.30902, -0.5, 0.30902, 0.80902, 0.30902, -0.80902, 0.5,
		0.80902, -0.5, -0.30902, 0.5, 0.30902, 0.80902, -0.30902, -0.80902, 0.5,
		0.5, 0.30902, -0.80902, 0.30902, 0.80902, 0.5, 0.80902, -0.5, 0.30902,
		-0.5, -0.30902, -0.80902, -0.30902, -0.80902, 0.5, -0.80902, 0.5, 0.30902,
		-0.30902, -0.80902, -0.5, 0.80902, -0.5, 0.30902, -0.5, -0.30902, 0.80902,
		0.30902, -0.80902, -0.5, 0.80902, 0.5, -0.30902, 0.5, -0.30902, 0.80902,
		-0.30902, 0.80902, -0.5, -0.80902, -0.5, -0.30902, -0.5, 0.30902, 0.80902,
		0.5, -0.30902, -0.80902, -0.30902, 0.80902, -0.5, 0.80902, 0.5, 0.30902,
		-0.5, 0.30902, -0.80902, 0.30902, -0.80902, -0.5, -0.80902, -0.5, 0.30902,
		-0.80902, -0.5, -0.30902, 0.5, -0.30902, -0.80902, 0.30902, -0.80902, 0.5,
		0.80902, 0.5, -0.30902, -0.5, 0.30902, -0.80902, -0.30902, 0.80902, 0.5,
		0.80902, -0.5, 0.30902, 0.5, 0.30902, -0.80902, 0.30902, 0.80902, 0.5,
		-0.80902, 0.5, 0.30902, -0.5, -0.30902, -0.80902, -0.30902, -0.80902, 0.5,
		-0.5, -0.30902, 0.80902, -0.30902, -0.80902, -0.5, 0.80902, -0.5, 0.30902,
		0.5, 0.30902, 0.80902, 0.30902, 0.80902, -0.5, -0.80902, 0.5, 0.30902,
		0, 0, 1, 1, 0, 0, 0, 1, 0,
		0, 1, 0, 0, 0, 1, 1, 0, 0,
		0, -1, 0, 0, 0, 1, -1, 0, 0,
		0, 0, -1, -1, 0, 0, 0, 1, 0,
		0, -1, 0, 0, 0, -1, 1, 0, 0,
		0, 1, 0, 0, 0, -1, -1, 0, 0,
		-0.80902, -0.5, 0.30902, -0.5, 0.30902, -0.80902, 0.30902, -0.80902, -0.5,
		0.80902, -0.5, -0.30902, -0.5, -0.30902, -0.80902, 0.30902, 0.80902, -0.5,
		0.5, 0.30902, -0.80902, -0.30902, -0.80902, -0.5, -0.80902, 0.5, -0.30902,
		-0.30902, -0.80902, -0.5, -0.80902, 0.5, -0.30902, 0.5, 0.30902, -0.80902,
		-0.80902, 0.5, -0.30902, 0.5, 0.30902, -0.80902, -0.30902, -0.80902, -0.5,
		-0.5, -0.30902, -0.80902, 0.30902, 0.80902, -0.5, 0.80902, -0.5, -0.30902,
		-0.5, 0.30902, -0.80902, -0.30902, 0.80902, 0.5, 0.80902, 0.5, -0.30902,
		0, 0, -1, 1, 0, 0, 0, -1, 0,
		-0.80902, 0.5, 0.30902, 0.5, 0.30902, 0.80902, 0.30902, 0.80902, -0.5,
		0.80902, 0.5, -0.30902, 0.5, -0.30902, 0.80902, 0.30902, -0.80902, -0.5,
		-0.30902, 0.80902, -0.5, 0.80902, 0.5, 0.30902, 0.5, -0.30902, -0.80902,
		0.5, -0.30902, -0.80902, 0.30902, -0.80902, 0.5, -0.80902, -0.5, -0.30902,
		-0.80902, -0.5, -0.30902, -0.5, 0.30902, 0.80902, -0.30902, 0.80902, -0.5,
		-0.30902, -0.80902, 0.5, -0.80902, 0.5, 0.30902, -0.5, -0.30902, -0.80902,
		-0.30902, 0.80902, 0.5, 0.80902, 0.5, -0.30902, -0.5, 0.30902, -0.80902,
		1, 0, 0, 0, -1, 0, 0, 0, -1,
		0.30902, 0.80902, -0.5, 0.80902, -0.5, -0.30902, -0.5, -0.30902, -0.80902,
		0.30902, -0.80902, -0.5, -0.80902, -0.5, 0.30902, -0.5, 0.30902, -0.80902,
		-1, 0, 0, 0, 1, 0, 0, 0, -1,
		0.80902, 0.5, 0.30902, 0.5, -0.30902, -0.80902, -0.30902, 0.80902, -0.5,
		0.30902, -0.80902, 0.5, -0.80902, -0.5, -0.30902, 0.5, -0.30902, -0.80902,
		-0.5, 0.30902, 0.80902, -0.30902, 0.80902, -0.5, -0.80902, -0.5, -0.30902,
		0, 0, 1, -1, 0, 0, 0, -1, 0,
		0.5, -0.30902, 0.80902, 0.30902, -0.80902, -0.5, 0.80902, 0.5, -0.30902,
		0.30902, 0.80902, 0.5, 0.80902, -0.5, 0.30902, 0.5, 0.30902, -0.80902,
		0.80902, -0.5, 0.30902, -0.5, -0.30902, 0.80902, -0.30902, -0.80902, -0.5,
		-0.5, -0.30902, 0.80902, 0.30902, 0.80902, 0.5, -0.80902, 0.5, -0.30902,
		0.5, 0.30902, 0.80902, -0.30902, -0.80902, 0.5, 0.80902, -0.5, -0.30902
	};

	int idx = n % 60;

	std::vector<float> matrix(12, 0);
	for (int r = 0; r < 3; ++r) {
		for (int c = 0; c < 3; ++c) {
			matrix[r*4 + c] = matrices[idx*9 + r*3 + c];
		}
	}

	Transform t3d(matrix);
	return t3d;
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
		if ( tmpaz > ( cap_sig/2.0f ) )tmpaz = cap_sig - tmpaz;

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



