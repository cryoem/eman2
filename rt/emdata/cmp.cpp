#include "emdata.h"
#include "util.h"
#include <stdio.h>

using namespace EMAN;

void mycmp(string cmpname, EMData& e1, EMData& e2)
{
	Dict params;
    float cmp1 = e1.cmp(cmpname, &e1, params);
    float cmp2 = e1.cmp(cmpname, &e1, params);
    
    printf("%12s: cmp(a, a) = %f; cmp(a, b) = %f\n", cmpname.c_str(), cmp1, cmp2);
}


int main(int argc, char* argv[])
{
    
    Util::set_log_level(argc, argv);
    
    char filename[1024];
    sprintf(filename, "%s/images/samesize1.mrc", getenv("HOME"));    

    EMData e1;
    e1.read_image(filename);

    sprintf(filename, "%s/images/samesize2.mrc", getenv("HOME"));    
    EMData e2;
    e2.read_image(filename);

    mycmp("Dot", e1, e2);
    mycmp("Variance", e1, e2);
    mycmp("Phase", e1, e2);
    mycmp("FRC", e1, e2);
    
    mycmp("HelloWorld", e1, e2);
    
    return 0;
}

	
