#include "emobject.h"
#include "emcache.h"
#include "emutil.h"

using namespace EMAN;

int test_emobject()
{
 
    EMObject e1 = EMObject();
    e1.get_int();

    EMObject e2 = EMObject(12.43);
    e2.get_int();
    e2.get_float();
    e2.get_string();
    
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
    
    assert(d["a"].get_int() == d.get("a").get_int());
    
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
	em->add(name, f);

	float* f1 = em->get(name);
	assert(f == f1);
    }

    return 0;
}
    

int main()
{
    int err = 0;

    err = test_emobject();
    err = test_emcache();
    err = test_Dict();
    
    return err;
}
