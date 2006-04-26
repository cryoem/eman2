/**
 * $Id$
 */
#include "emdata.h"
#include "processor.h"
#include "cmp.h"
#include "aligner.h"
#include "projector.h"

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


EMData* EMData::process(const string & processorname, const Dict & params)
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

