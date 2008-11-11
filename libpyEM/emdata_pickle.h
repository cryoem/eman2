/*
 * Author: Grant Tang, 07/12/2006 (sludtke@bcm.edu)
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

#ifndef EMDATA_PICKLE_H_
#define EMDATA_PICKLE_H_

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/back_reference.hpp>

#include "emdata.h"
#include "emobject.h"
#include "typeconverter.h"

//using namespace boost;

struct EMData_pickle_suite : boost::python::pickle_suite
{
/*	static
	boost::python::tuple
	getinitargs(const EMAN::EMData& em)
	{
		using namespace boost::python;
		return make_tuple();
	}
*/
	static
	boost::python::tuple
	getstate(boost::python::object em_obj)
	{
		using namespace boost::python;
		EMAN::EMData const& em = extract<EMAN::EMData const&>(em_obj)();
		
		return make_tuple(em_obj.attr("__dict__"),
							em.get_flags(), 
							em.get_changecount(),
							em.get_xsize(), 
							em.get_ysize(), 
							em.get_zsize(),
							em.get_xoff(),
							em.get_yoff(),
							em.get_zoff(),
							em.get_path(),
							em.get_pathnum(),
							em.get_attr_dict(),
							em.get_translation(),
							em.get_data_pickle(),
							em.get_supp_pickle());
	}
	
	static
	void
	setstate(boost::python::object em_obj, boost::python::tuple state)
	{
		using namespace boost::python;
		EMAN::EMData & em = extract<EMAN::EMData &>(em_obj)();
		
		if(len(state) != 15) {
			PyErr_SetObject( PyExc_ValueError, 
							("expected 15-item tuple in call to__set_state__; got %s"
							 % state).ptr()
							 );
			throw_error_already_set();
		}
		
		dict d = extract<dict>(em_obj.attr("__dict__"))();
		d.update(state[0]);
		
		int flags = extract<int>(state[1]);
		em.set_flags(flags);
		
		int changecount = extract<int>(state[2]);
		em.set_changecount(changecount);
		
		int nx = extract<int>(state[3]);
		int ny = extract<int>(state[4]);
		int nz = extract<int>(state[5]);
		em.set_size(nx, ny, nz);
		
		int xoff = extract<int>(state[6]);
		int yoff = extract<int>(state[7]);
		int zoff = extract<int>(state[8]);
		em.set_xyzoff(xoff, yoff, zoff);
		
		string path = extract<string>(state[9]);
		em.set_path(path);
		
		int pathnum = extract<int>(state[10]);
		em.set_pathnum(pathnum);
		
		EMAN::Dict attr_dict = extract<EMAN::Dict>(state[11]);
		em.set_attr_dict(attr_dict);
		
		EMAN::Vec3f all_translation = extract<EMAN::Vec3f>(state[12]);
		em.set_translation(all_translation);
		
		//vector<float> vf = extract< vector<float> >(state[13]);
		std::string vf = extract< std::string >(state[13]);
		em.set_data_pickle(vf);
		
		int fake_supp = extract<int>(state[14]);
		em.set_supp_pickle(fake_supp);
	}
	
	static bool getstate_manage_dict() {return true;}
};


struct Vec3f_pickle_suite : boost::python::pickle_suite
{
	static
	boost::python::tuple
	getinitargs(const EMAN::Vec3f& v)
	{
		using namespace boost::python;
		return make_tuple(v[0], v[1], v[2]);
	}
};

struct Transform_pickle_suite : boost::python::pickle_suite
{
	static
	boost::python::tuple
	getinitargs(const EMAN::Transform& t)
	{
		using namespace boost::python;
		return make_tuple(t.get_matrix());
	}
};

struct EMAN1Ctf_pickle_suite : boost::python::pickle_suite
{
	static
		boost::python::tuple
		getinitargs(const EMAN::EMAN1Ctf& ctf)
		{
			using namespace boost::python;
			return make_tuple(ctf.to_vector());
		}
};

struct EMAN2Ctf_pickle_suite : boost::python::pickle_suite
{
	static
		boost::python::tuple
		getinitargs(const EMAN::EMAN2Ctf& ctf)
		{
			using namespace boost::python;
			return make_tuple(ctf.to_vector());
		}
};

#endif	//EMDATA_PICKLE_H_
