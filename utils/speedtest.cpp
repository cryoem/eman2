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

/* This program runs a set of speed tests on the current machine

 * This program runs a set of speed tests in the current computer. It should
 * have at least 1 unloaded processor, and minimal i/o occuring when this
 *  is run. It will give a single value which should be generally proportional
 * to how fast refine will run on a given computer. It is a 'real world' test,
 * in that it simulates a variety of actual operations performed by a
 * refinement procedure. Don't compare values given by speed test in
 * different versions of EMAN, since the underlying routines may be
 * improved with time.
 *
 * */

#include <ctime>

#include "emdata.h"
#include <sstream>
using std::stringstream;

#include <iostream>
using std::cout;
using std::endl;

using namespace EMAN;

#define CPS	(CLOCKS_PER_SEC)


int main(int argc, char *argv[])
{
	int SIZE = 96;
    int NTT = 500;

    int slow = 0;
    int low = 0;
    int big = 0;
    int newali = 0;
    int vg = 0;
#ifdef EMAN2_USING_CUDA
    int gpu = 0;
#endif

    if (argc > 1) {
		if (Util::sstrncmp(argv[1], "slow"))
			slow = 1;
		else if (Util::sstrncmp(argv[1], "best"))
			slow = 3;
		else if (Util::sstrncmp(argv[1], "low"))
			low = 1;
		else if (Util::sstrncmp(argv[1], "refine"))
			newali = 1;
		else if (Util::sstrncmp(argv[1], "big"))
			big = 1;
		else if (Util::sstrncmp(argv[1], "valgrind"))
			vg=1;
#ifdef EMAN2_USING_CUDA
		else if (Util::sstrncmp(argv[1], "gpu")) {
			cout << "GPU mode selected" << endl;
			gpu=1;
		}
#endif
    }
    if (argc > 2) {
    	// assume it's the boxsize
    	stringstream ss;
    	ss << argv[2];
    	ss >> SIZE;
    	if (SIZE < 0 || SIZE > 1024) {
    		cout << "The boxsize " << SIZE << " is invalid. Please specify a box size in the range of (0,1024]" << endl;
    		throw;
    	}
    }

    // This is some testing code used to run valgrind tests on
    // parts of the core library. It is not for benchmarking
    // purposes.
    if (vg) {
		EMData *a = new EMData();
		a->set_size(128,128,1);
		a->process_inplace("testimage.scurve");
		EMData *b = a->rot_scale_trans2D(0.0,1.0,1.0,0.0);
		EMData *d = b->align("rotate_translate",a,Dict(),"sqeuclidean");
		Dict r=b->get_attr_dict();
		printf("%p\n",&r);
		printf("nx = %d\n",(int)r["nx"]);
		delete a;
		delete b;
		delete d;
		exit(0);
    }

    printf("This could take a few minutes. Please be patient.\n");
    if (big) {
		NTT = 100;
		SIZE = 320;
    }

    printf("Initializing\n");

    EMData pat;
    pat.set_size(SIZE, SIZE, 1);
    float *d = pat.get_data();

    for (int i = -SIZE / 2; i < SIZE / 2; i++) {
		for (int j = -SIZE / 2; j < SIZE / 2; j++) {
			d[(i + SIZE / 2) + (j + SIZE / 2) * SIZE] =
				-3.0f * exp(-Util::square((fabs((float)i) + fabs((float)j))) / 10.0f) +
				exp(-Util::square((fabs((float)i) + fabs((float)j) / 2.0f)) / 100.0f) *
				(abs(i) < 2 ? 2.0f : 1.0f);
		}
    }
    pat.update();
    pat.process_inplace("normalize.circlemean");
    pat.process_inplace("mask.sharp", Dict("outer_radius", pat.get_xsize()/2));

    EMData *data[5000];

    for (int i = 0; i < NTT; i++) {
		data[i] = pat.copy();
		float *d = data[i]->get_data();

		for (int j = 0; j < SIZE * SIZE; j++) {
			d[j] += Util::get_gauss_rand(0, 1.0);
		}
		data[i]->update();
		data[i]->process_inplace("normalize.circlemean");
		data[i]->process_inplace("mask.sharp", Dict("outer_radius", data[i]->get_xsize()/2));

		if (i < 5) {
			data[i]->write_image("speed.hed", i, EMUtil::IMAGE_IMAGIC);
		}
    }

    if (low) {
		printf("Low level tests starting. Please note that compiling with optimization may invalidate certain tests. Also note that overhead is difficult to compensate for, so in most cases it is not dealt with.\n");

		float t1 = (float) clock();
		for (float fj = 0; fj < 500.0; fj += 1.0) {
			for (int i = 0; i < NTT / 2.0; i += 1)
				data[i]->dot(data[i + NTT / 2]);
		}

		float t2 = (float) clock();
		float ti = (t2 - t1) / (float) CPS;
		printf("Baseline 1: %d, %d x %d dot()s in %1.1f sec  %1.1f/sec-> ~%1.2f mflops\n",
			   500 * NTT / 2, SIZE, SIZE, ti, 500 * NTT / (2.0 * ti),
			   SIZE * SIZE * 4 * 500.0 * NTT / (1000000.0 * ti));

		t1 = (float) clock();
		for (float fj = 0; fj < 500.0; fj += 1.0) {
			for (float fi = 0; fi < NTT / 2.0; fi += 1.0)
				data[1]->dot(data[12]);
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 1a: %d, %d x %d optimized cached dot()s in %1.1f s %1.1f/s-> ~%1.2f mflops\n",
			   500 * NTT / 2, SIZE, SIZE, ti, 500 * NTT / (2.0 * ti),
			   SIZE * SIZE * 4 * 500.0 * NTT / (1000000.0 * ti));

		t1 = (float) clock();
		for (int j = 0; j < 500; j++) {
			for (int i = 0; i < NTT / 2; i++) {
				Dict d;
				d["keepzero"] = 1;
				data[i]->cmp("sqeuclidean", data[i + NTT / 2], d);
			}
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 2: %d, %d x %d lcmp()s in %1.1f sec   %1.1f/sec\n", 500 * NTT / 2, SIZE,
			   SIZE, ti, 500 * NTT / (2.0 * ti));

		t1 = (float) clock();
		for (int j = 0; j < 100; j++) {
			for (int i = 0; i < NTT / 2; i++) {
				Dict params;
				data[i]->cmp("phase", data[i + NTT / 2], params);
			}
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 3: %d, %d x %d pcmp()s in %1.1f sec  %1.1f/sec\n", 100 * NTT / 2, SIZE,
			   SIZE, ti, 100 * NTT / (2.0 * ti));

		t1 = (float) clock();
		for (int j = 0; j < 100; j++) {
			for (int i = 0; i < NTT / 2; i++) {
				Dict params;
				data[i]->cmp("frc", data[i + NTT / 2], params);
			}
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 4: %d, %d x %d fscmp()s in %1.1f sec  %1.1f/sec\n", 100 * NTT / 2, SIZE,
			   SIZE, ti, 100 * NTT / (2.0 * ti));

		t1 = (float) clock();
		for (int j = 0; j < 500; j++) {
			for (int i = 0; i < NTT / 2; i++)
				data[i]->process_inplace("math.absvalue");
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 5a: %d, %d x %d fabs in %1.1f sec -> ~%1.2f mfabs/sec\n", 500 * NTT / 2,
			   SIZE, SIZE, ti, SIZE * SIZE * 500.0 * NTT / (1000000.0 * ti));

		t1 = (float) clock();
		for (int j = 0; j < 100; j++) {
			for (int i = 0; i < NTT / 2; i++)
				data[i]->process_inplace("math.sqrt");
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 5b: %d, %d x %d sqrts in %1.1f sec -> ~%1.2f msqrt/sec\n", 100 * NTT / 2,
			   SIZE, SIZE, ti, SIZE * SIZE * 100.0 * NTT / (1000000.0 * ti));

		d = data[0]->get_data();
		t1 = (float) clock();
		for (int j = 0; j < 100; j++) {
			for (int i = 0; i < NTT / 2; i++) {
				for (int k = 0; k < SIZE * SIZE; k++)
					(void)sqrt(d[k]);
			}
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 5c: %d, %d x %d sqrts in %1.1f sec -> ~%1.2f msqrt/sec (cached)\n",
			   100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 100.0 * NTT / (1000000.0 * ti));
		data[0]->update();

		d = data[0]->get_data();
		t1 = (float) clock();
		for (int j = 0; j < 100; j++) {
			for (int i = 0; i < NTT / 2; i++) {
				for (int k = 0; k < SIZE * SIZE; k++)
					(void)cos(d[k]);
			}
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 5d: %d, %d x %d cos in %1.1f sec -> ~%1.2f mcos/sec (cached)\n",
			   100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 100.0 * NTT / (1000000.0 * ti));
		data[0]->update();

		d = data[0]->get_data();
		t1 = (float) clock();
		for (int j = 0; j < 100; j++) {
			for (int i = 0; i < NTT / 2; i++) {
				for (int k = 0; k < SIZE * SIZE - 1; k++)
					(void)hypot(d[k], d[k + 1]);
			}
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 5e: %d, %d x %d hypot in %1.1f sec -> ~%1.2f mhypot/sec (cached)\n",
			   100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 100.0 * NTT / (1000000.0 * ti));
		data[0]->update();

		d = data[0]->get_data();
		t1 = (float) clock();
		for (int j = 0; j < 1000; j++) {
			for (int i = 0; i < NTT / 2; i++) {
				for (int k = 0; k < SIZE * SIZE - 1; k++) {
					float f = d[k] * d[k + 1];
					f = f + f;
				}
			}
		}

		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 5f: %d, %d x %d mult in %1.1f sec -> ~%1.2f mmult/sec (cached)\n",
			   100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 1000.0 * NTT / (1000000.0 * ti));
		data[0]->update();

		d = data[0]->get_data();
		t1 = (float) clock();
		for (int j = 0; j < 500; j++) {
			for (int i = 0; i < NTT / 2; i++) {
				for (int k = 0; k < SIZE * SIZE - 1; k++) {
					float a = d[k] / d[k + 1];
					a = a + a;
				}
			}
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 5g: %d, %d x %d div in %1.1f sec -> ~%1.2f mdiv/sec (cached)\n",
			   100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 500.0 * NTT / (1000000.0 * ti));
		data[0]->update();

		d = data[0]->get_data();
		t1 = (float) clock();
		for (int j = 0; j < 500; j++) {
			for (int i = 0; i < NTT / 2; i++) {
				for (int k = 0; k < SIZE * SIZE - 1; k++) {
					float f = fabs(d[k]);
					f = f + f;
				}
			}
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 5h: %d, %d x %d fabs in %1.1f sec -> ~%1.2f fabs/sec (cached)\n",
			   100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 500.0 * NTT / (1000000.0 * ti));
		data[0]->update();

		d = data[0]->get_data();
		t1 = (float) clock();
		for (int j = 0; j < 500; j++) {
			for (int i = 0; i < NTT / 2; i++) {
				for (int k = 0; k < SIZE * SIZE - 1; k++) {
					(void)atan2(d[k], d[k + 1]);
					(void)hypot(d[k], d[k + 1]);
				}
			}
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 5i: %d, %d x %d ri2ap in %1.1f sec -> ~%1.2f ri2ap/sec (cached)\n",
			   100 * NTT / 2, SIZE, SIZE, ti, SIZE * SIZE * 500.0 * NTT / (1000000.0 * ti));

		data[0]->update();
		t1 = (float) clock();

		for (int i = 0; i < NTT * 100; i++) {
			EMData *cp = data[i % NTT]->copy();
			Dict d("n",2);
			cp->process_inplace("math.meanshrink",d);
			if( cp )
			{
				delete cp;
				cp = 0;
			}
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 6:  %1.1f sec %f meanshrink x 2/sec\n", ti, NTT * 100.0 / ti);

		EMData *d1a = data[0]->copy();
		t1 = (float) clock();

		for (int i = 0; i < NTT * 1000; i++) {
			d1a->do_fft();
			d1a->update();
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 7:  %1.1f sec %f ffts/sec\n", ti, NTT * 1000 / ti);

		d1a = d1a->copy();
		t1 = (float) clock();

		for (int i = 0; i < NTT * 1000; i++) {
			d1a->translate(-1, -3, 0);
		}
		t2 = (float) clock();
		ti = (t2 - t1) / (float) CPS;
		printf("Baseline 8:  %1.1f sec   %f translates/sec\n", ti, NTT * 1000 / ti);

		return 0;
    }


    EMData *tmp = 0;
    float t1 = (float) clock();
    for (int i = 0; i < 3; i++) {
#ifdef EMAN2_USING_CUDA
		if (gpu) {
// 			cout << "Using gpu" << endl;
			data[i]->set_gpu_rw_current();
			data[i]->cuda_lock();
		}
#endif
		for (int j = 5; j < NTT; j++) {
			if (slow == 3) {
				tmp = data[i]->align("rtf_best", data[j],
									 Dict("flip", (EMData*)0,"maxshift", SIZE/8));
			}
			else if (slow == 1) {
				Dict d;
				d["flip"] = (EMData*)0;
				d["maxshift"] = SIZE/8;
				tmp = data[i]->align("rtf_slow", data[j], d);
			}
			else if (newali == 1) {
				tmp = data[i]->align("rotate_translate_flip", data[j], Dict());
				EMData* tmp2 = tmp->align("refine", data[j], Dict());
				delete tmp2;
				tmp2 = 0;
			}
			else {
#ifdef EMAN2_USING_CUDA
				if (gpu) {
					data[j]->set_gpu_rw_current();
				}
#endif
				tmp = data[i]->align("rotate_translate_flip", data[j], Dict());
			}
	    	if( tmp )
	    	{
				delete tmp;
				tmp = 0;
	    	}
			if (j % 10 == 0) {
				putchar('.');
				fflush(stdout);
			}
		}
#ifdef EMAN2_USING_CUDA
		if (gpu) {
			data[i]->cuda_unlock();
		}
#endif
		putchar('\n');
    }
    float t2 = (float) clock();

    float ti = (t2 - t1) / (float) CPS;
    if (slow == 2)
		ti *= 10.0;
    if (!big && !slow && !newali) {
		printf("\nFor comparison (note these numbers may change from release to release)\n");
		printf("An AMD Athlon (32 bit) 900Mhz SF --------------------------------  360\n");
		printf("An AMD Athlon XP 2400+ (32 bit) 2.0Ghz SF ----------------------- 1010\n");
		printf("An AMD Athlon XP 2600+ (32 bit) 2.0Ghz SF ----------------------- 1090\n");
		printf("An AMD Athlon 64 3700+ 2.2Ghz SF -------------------------------- 1530\n");
		printf("An AMD Athlon 64 X2 3800+ 2.0Ghz SF ----------------------------- 1578\n");
		printf("An AMD Athlon 64 FX-51 2.2Ghz SF -------------------------------- 1760\n");
		printf("An AMD Opteron 248 2.2Ghz SF ------------------------------------ 1870\n");
		printf("An Intel Core2 T7200 2.0Ghz SF ---------------------------------- 1990\n");
		printf("An Intel Xeon E5335 2.0Ghz SF ----------------------------------- 2010\n");
		printf("An AMD Opteron 280 2.4Ghz SF ------------------------------------ 2130\n");
		printf("An Intel Core2 6700 2.66Ghz SF ---------------------------------- 2600\n");
		printf("An Intel Core2 Duo T9400 2.53Ghz SF ----------------------------- 2730\n");
		printf("An Intel Xeon E5430 2.66Ghz SF ---------------------------------- 2800\n");
		printf("An Intel Xeon X5355 2.66Ghz SF ---------------------------------- 2920\n");
		printf("An Intel Xeon X5450 3.0Ghz SF ----------------------------------- 3220\n");
		printf("An Intel Xeon X5460 3.16Ghz SF ---------------------------------- 3320\n");
		printf("\nYour machines speed factor = %1.1f\n", 25000.0 / ti);
		printf("\n\nThis repesents %1.2f RTFAligns/sec\n",
			   3.0 * ((slow == 2 ? NTT / 10 : NTT) - 5) / ti);
    }
    else if (big && !slow) {
		printf("\nYour machines speed factor = %1.1f\n", 72000.0 / ti);
    }
    else {
		printf("\nYour machines speed factor on this test = %1.1f\n", 25000.0 / ti);
    }
    return 0;
}


