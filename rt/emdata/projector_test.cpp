// test projector

#include "ctf.h"
#include "log.h"
#include "emdata.h"
#include "xydata.h"
#include "assert.h"
#include "projector.h"

#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <sys/time.h>

#ifndef CLK_TCK
#define CLK_TCK 60
#endif

double mytimer()
{
	struct tms use;
	double tmp;
	times(&use);
	tmp = use.tms_utime;
	tmp += use.tms_stime;
	return (double)(tmp) / CLK_TCK;
}

using namespace EMAN;

using std::cout;
using std::cin;
using std::endl;
using std::string;

//#include "spidutil.h"

int main()
{
	EMData *volume = new EMData(); // initial volume
	Dict myparams;
	
	vector<float> angles;
	float delta, tlb, tub, plb, pub;
	int ndim, nx, ny, nz, nang, status;
	int nref;
	double t0, t1;
	
	cout << "testing projector " << endl;
	volume->read_image("vol001.tfc");
	
	ndim = volume->get_ndim();
	nx = volume->get_xsize();
	ny = volume->get_ysize();
	nz = volume->get_zsize();

	// print some stats
	printf("ndim = %d\n", ndim);
	printf("nx = %d\n", nx);
	printf("ny = %d\n", ny);
	printf("nz = %d\n", nz);
	
	float mean = volume->get_attr("mean");
	cout << " mean = " << mean << endl;
	
	float min = volume->get_attr("minimum");
	cout << " min = " << min << endl;
	
	float max = volume->get_attr("maximum");
	cout << " max = " << max << endl;

	// do a single projection
	myparams["angletype"] = "SPIDER";
	myparams["phi"] = 10.0;
	myparams["theta"] = 30.0;
	myparams["psi"] = 50.0;
	
	EMData* proj = volume->project("pawel", myparams);

	mean = proj->get_attr("mean");
	cout << " projmean = " << mean << endl;

	proj->write_image("proj1.spi",0,EMUtil::IMAGE_SINGLE_SPIDER);

	delta = 15.0;
	tlb = 0.0;
	tub = 90.0;
	plb = 0.0;
	pub = 359.9;

	angles = Util::voea(delta, tlb, tub, plb, pub);
	nang = angles.size()/3;
	printf(" size(angles) = %d\n", nang);

	myparams["anglelist"] = angles;
	myparams["angletype"] = "SPIDER";
	myparams["radius"] = 30.0;
	t0 = mytimer();
	EMData* refprj = volume->project("pawel", myparams);
	t1 = mytimer() - t0;
	nref = refprj->get_zsize();
	printf(" nref = %d\n", nref);
	
	printf("Projection time = %11.3e\n", t1);
	// refprj->write_image("refprj.img", -1);

	status = 0;

	if (volume)	delete volume;
	if (proj)	delete proj;
	if (refprj)	delete refprj;
	return status;
}
