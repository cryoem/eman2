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

#include "emobject.h"
#include "emcache.h"
#include "emutil.h"
#include "emdata.h"
#include "xydata.h"
#include "emassert.h"
#include "processor.h"

using namespace EMAN;

int test_emobject()
{
 
    EMObject e1 = EMObject();
    int n = e1;
    
    EMObject e2 = EMObject(12.43);
    n = e2;
    float f = e2;
    f = 2;
    const char *s1 = e2;
    Assert(strcmp(s1, "") == 0);

    EMData * image1 = new EMData();
    EMObject e3(image1);
    EMData * image2 = e3;
    Assert(image1 == image2);

    XYData * xy1 = new XYData();
    EMObject e4(xy1);
    XYData * xy2 = e4;
    Assert(xy1 == xy2);
    
    if( image1 )
    {
    	delete image1;
   	 	image1 = 0;
    }
    
    return 0;
}

int test_Dict()
{
    Dict d;
    d["a"] = 1;
    d["b"] = 2;
    d["c"] = 3;

    bool f1 = d.has_key("hello");
    bool f2 = d.has_key("c");

    Assert(f1 == false);
    Assert(f2 == true);
    
    Assert((int)d["a"] == (int)d.get("a"));
    
    EMUtil::dump_dict(d);

    return 0;
}


int test_emcache()
{
    EMCache<float>* em = new EMCache<float>(4);

    for (int i = 0; i < 20; i++) {
	char name[32];
	sprintf(name, "name%d", i);

	float* f = new float(i);
	em->add(string(name), f);

	float* f1 = em->get(name);
	Assert(f == f1);
    }

    return 0;
}


void test_factory()
{
    //Factory<Processor>::instance()->dump();
    EMAN::dump_processors();
}


int main()
{
    int err = 0;

    err = test_emobject();
    err = test_emcache();
    err = test_Dict();

    test_factory();
    
    return err;
}
