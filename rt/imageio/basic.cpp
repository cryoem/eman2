#include "emobject.h"
#include "emcache.h"
#include "emutil.h"
#include "emdata.h"
#include "xydata.h"
#include <assert.h>
#include "filter.h"

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
    assert(strcmp(s1, "") == 0);

    EMData * image1 = new EMData();
    EMObject e3(image1);
    EMData * image2 = e3;
    assert(image1 == image2);

    XYData * xy1 = new XYData();
    EMObject e4(xy1);
    XYData * xy2 = e4;
    assert(xy1 == xy2);
    
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

    assert(f1 == false);
    assert(f2 == true);
    
    assert((int)d["a"] == (int)d.get("a"));
    
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
	assert(f == f1);
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
