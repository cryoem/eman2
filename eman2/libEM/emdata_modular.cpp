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

#include "emdata.h"
#include "processor.h"
#include "cmp.h"
#include "aligner.h"
#include "projector.h"
#include "analyzer.h"

using namespace EMAN;

void EMData::process_inplace(const string & processorname, const Dict & params)
{
	ENTERFUNC;
	Processor *f = Factory < Processor >::get(processorname, params);
	if (f) {
		f->process_inplace(this);
		if( f )
		{
			delete f;
			f = 0;
		}
	}
	EXITFUNC;
}

void EMData::process_inplace(Processor * p)
{
	ENTERFUNC;
	if(p) {
		p->process_inplace(this);
	}
	EXITFUNC;
}

EMData* EMData::process(const string & processorname, const Dict & params) const
{
	ENTERFUNC;
	Processor *f = Factory < Processor >::get(processorname, params);
	EMData * result = 0;
	if (f) {
		result = f->process(this);
		if( f )
		{
			delete f;
			f = 0;
		}
	}
	return result;
	EXITFUNC;
}

EMData * EMData::process(Processor * p) const
{
	ENTERFUNC;
	EMData * result = 0;
	if(p) {
		result = p->process(this);
	}
	return result;
	EXITFUNC;
}

float EMData::cmp(const string & cmpname, EMData * with, const Dict & params)
{
	ENTERFUNC;
	float result = 0;
	Cmp *c = Factory < Cmp >::get(cmpname, params);
	if (c) {
		result = c->cmp(this, with);
		if( c )
		{
			delete c;
			c = 0;
		}
	}

	EXITFUNC;
	return result;
}


EMData *EMData::align(const string & aligner_name, EMData * to_img,
					  const Dict & params, const string & cmp_name, const Dict& cmp_params)
{
	ENTERFUNC;
	EMData *result = 0;
	Aligner *a = Factory < Aligner >::get(aligner_name, params);
	if (a) {
		if (cmp_name == "") {
			result = a->align(this, to_img);
		}
		else {
			result = a->align(this, to_img, cmp_name, cmp_params);
		}
		if( a )
		{
			delete a;
			a = 0;
		}
	}

	EXITFUNC;
	return result;
}

vector<Dict> EMData::xform_align_nbest(const string & aligner_name, EMData * to_img,
								const Dict & params, const unsigned int nsoln, const string & cmp_name,
										const Dict& cmp_params)
{
	ENTERFUNC;
	Aligner *a = Factory < Aligner >::get(aligner_name, params);
	vector<Dict> result;
	if (a) {
		result = a->xform_align_nbest(this,to_img,nsoln,cmp_name,cmp_params);
	}

	return result;
}

EMData *EMData::project(const string & projector_name, const Dict & params)
{
	ENTERFUNC;
	EMData *result = 0;
	Projector *p = Factory < Projector >::get(projector_name, params);
	if (p) {
		result = p->project3d(this);
		if( p )
		{
			delete p;
			p = 0;
		}
	}

	EXITFUNC;
	return result;
}


EMData *EMData::project(const string & projector_name, const Transform & t3d)
{
	ENTERFUNC;
	EMData *result = 0;
	Dict params;
	params["transform"] = (Transform*) &t3d;
	Projector *p = Factory < Projector >::get(projector_name, params);
	if (p) {
		result = p->project3d(this);
		if( p )
		{
			delete p;
			p = 0;
		}
	}

	EXITFUNC;
	return result;
}

EMData *EMData::backproject(const string & projector_name, const Dict & params)
{
	ENTERFUNC;
	EMData *result = 0;
	Projector *p = Factory < Projector >::get(projector_name, params);
	if (p) {
		result = p->backproject3d(this);
		if( p )
		{
			delete p;
			p = 0;
		}
	}

	EXITFUNC;
	return result;
}
