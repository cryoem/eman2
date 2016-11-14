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
#include "ctf.h"

#include <cctype>

using namespace EMAN;
using std::map;

void usage(const char *progname)
{
    printf("\n%s inputfile outputfile [-vN] [first=<n>] [last=<n>] [inplace] [edgenorm] [norm[=<mean>,<sigma>]] [invert] [flip] [rot=<angle>] [trans=<dx>,<dy>] [center] [acfcenter] [scale=<sca>] [clip=<x,y>] [shrink=<n>] [meanshrink=<n>] [apix=<A/pix>] [lp=<filt r>] [hp=<filt r>]  [blfilt=<sigma1,sigma2,iter,localwidth>] [mask=<radius>] [imask=<radius>] [noisemask]  [mrc]  [automask=x] [pif] [hdf] [png] [em] [spider] [spider-single] [pgm[=<low>,<hi>]] [vtk] [sharphp=<pixels>] [norefs] [rotav] [average]  [split=n] [ctfsplit] [interlv=<file 2 interleave>] [sym=<Cn>] [plt[=<plt file>]] [setsfpairs] [addnoise=<level>] [randomize=<dr>,<da>,<flip>]  [selfcl[=<steps>][,<mode>]] [radon] [rfp] [mraprep] [phot] [rsub] rfilt=<filtername><:param1=value1><...>\n", progname);
}

int main(int argc, char *argv[])
{

    if (argc < 3) {
		usage(argv[0]);
		exit(0);
    }

    Util::set_log_level(argc, argv);

    const int MAXMICROCTF = 1000;
    float *defocus_val = new float[MAXMICROCTF];
    float *bfactor_val = new float[MAXMICROCTF];

    string center = "center";
    string acfcenter = "acfcenter";
    string invert = "invert";
    string phot = "phot";

    string avgonly = "average";
    EMData *average = 0;

    string rfp = "rfp";
    string noisemask = "noisemask";
    string mraprep = "mraprep";
    string ctfsplit = "ctfsplit";
    string setsfpairs = "setsfpairs";
    string norefs = "norefs";
    string mrc8bit = "mrc8bit";

    string mrc = "mrc";
    string pif = "pif";
    string png = "png";
    string hdf = "hdf";
	string vtk = "vtk";
	
    string em = "em";
    string spidersingle = "spider-single";

    string spider = "spider";
    string rotav = "rotav";
    string inplace = "inplace";
    string rsub = "rsub";

    string flip = "flip";
    string edgenorm = "edgenorm";
    string radon = "radon";

    float blsigma1 = 0;
    float blsigma2 = 0;
    int bliter = 0;
    int blwidth = 0;

    map < string, bool > argdict;
    map < string, string > argfilters;

    float apix = 0;
    float lp = 0;
    float tlp = 0;

    float hp = 0;
    float sharphp = 0;

    int norm = 0;
    float nmean = 0;
    float nsig = 0;

    float ramp = 0;
    int mask = 0;
    int imask = 0;
    float automask = 0;
    
    float anoise = 0;
    float rot = 0;
    float dx = 0;
    float dy = 0;

    int rize = 0;
    int rizef = 0;
    float rizedx = 0;
    float rizeda = 0;

    float scale = 1;
    int clipx = -1;
    int clipy = -1;

    int shrink = 1;
    int csym = 0;

    int scl = 0;
    int sclmd = 0;

    int pgm = 0;
    int pgmlo = 0;
    int pgmhi = 0;

    bool rfilt = false;
    
    const char *interleaved_file = 0;

    int n0 = 0;
    int n1 = -1;

    argfilters[noisemask] = "MakeRadiusSquared";
    string filtername;
    Dict params_dict;
    
    for (int i = 3; i < argc; i++) {
		const char *arg = argv[i];

		if (Util::get_str_float(arg, "apix=", &apix)) {
		}
		else if (Util::get_str_float(arg, "lp=", &lp)) {
		}
		else if (Util::get_str_float(arg, "tlp=", &tlp)) {
		}
		else if (Util::get_str_float(arg, "hp=", &hp)) {
		}
		else if (Util::get_str_float(arg, "sharphp=", &sharphp)) {
		}
		else if (Util::get_str_int(arg, "first=", &n0)) {
		}
		else if (Util::get_str_int(arg, "last=", &n1)) {
		}
		else if (Util::get_str_float(arg, "ramp=", &ramp)) {
		}
		else if (Util::get_str_int(arg, "mask=", &mask)) {
		}
		else if (Util::get_str_int(arg, "imask=", &imask)) {
		}
		else if (Util::get_str_float(arg, "addnoise=", &anoise)) {
		}
		else if (Util::get_str_float(arg, "rot=", &rot)) {
		}
		else if (Util::get_str_float(arg, "automask=", &automask)) {
		}
		else if (Util::get_str_float(arg, "trans=", &dx, &dy)) {
		}
		else if (Util::get_str_float(arg, "scale=", &scale)) {
		}
		else if (Util::get_str_int(arg, "clip=", &clipx, &clipy)) {
		}
		else if (Util::get_str_int(arg, "shrink=", &shrink)) {
		}
		else if (strncmp(argv[i], "blfilt=", 7) == 0) {
			sscanf(argv[i] + 7, "%f,%f,%d,%d", &blsigma1, &blsigma2, &bliter, &blwidth);
		}
		else if (strncmp(arg, "sym=", 4) == 0) {
			if (tolower(argv[i][4]) != 'c') {
				LOGERR("Only C symmetry currently supported!\n");
				exit(1);
			}
			csym = atoi(&argv[i][5]);
		}
		else if (strncmp(arg, "randomize=", 10) == 0) {
			rize = 1;
			sscanf(&argv[i][10], "%f,%f,%d", &rizedx, &rizeda, &rizef);
			rizeda *= M_PI / 180.0f;
		}
		else if (Util::get_str_float(arg, "norm", &norm, &nmean, &nsig)) {
			if (norm != 2) {
				norm = 1;
			}
		}
		else if (strncmp(argv[i], "selfcl", 6) == 0) {
			if (argv[i][6] == '=') {
				sscanf(argv[i] + 7, "%d,%d", &scl, &sclmd);
				scl /= 2;
			}
			else {
				scl = 90;
			}
		}
		else if (strncmp(arg, "interlv=", 8) == 0) {
			interleaved_file = arg + 8;
		}
		else if (Util::get_str_int(arg, "pgm", &pgm, &pgmlo, &pgmhi)) {
			if (pgm != 2) {
				pgm = 1;
			}
		}
		else if (Util::sstrncmp(arg, noisemask.c_str())) {
			argfilters[noisemask] = "NoiseMask";
		}
		else if (strncmp(arg, "rfilt=", 6) == 0) {
			rfilt = true;
			string s(&arg[6]);
			size_t cpos = s.find(":");
	    
			if (cpos == string::npos) {
				filtername = s;
			}
			else {
				filtername = s.substr(0, cpos);
				size_t slen = s.length();
		
				while (cpos < slen) {
					size_t eqpos = s.find("=", cpos);
					string param1 = s.substr(cpos+1, eqpos-cpos-1);
		    
					cpos = s.find(":", cpos+1);
					string val1;
					if (cpos == string::npos) {
						val1 = s.substr(eqpos+1);
						cpos = slen;
					}
					else {
						val1 = s.substr(eqpos+1, cpos-eqpos-1);
					}
		    
					params_dict[param1] = atof(val1.c_str());
				}
			}
	    
		}
		else {
			argdict[arg] = true;
		}
    }



    EMData *d = new EMData();

    int nimg = EMUtil::get_image_count(argv[1]);
    if (nimg <= n1 || n1 < 0) {
		n1 = nimg - 1;
    }
    
    char outfile[1024];
    sprintf(outfile, "%s", argv[2]);
    char spiderformat[256];

    EMData *ld = new EMData();
    vector < float >sfcurve1;
    try {
		for (int i = n0; i <= n1; i++) {
			d->read_image(argv[1], i);
			int nx = d->get_xsize();
	
			if ((argdict[ctfsplit] != 0) && (i == n0 || !EMUtil::is_same_ctf(d, ld))) {
				Ctf *ctf = d->get_ctf();

				int j = 0;
				for (j = 1; j < argdict[ctfsplit]; j++) {
					if (defocus_val[j] == ctf->get_defocus() && bfactor_val[j] == ctf->get_bfactor()) {
						break;
					}
				}
				if (argdict[ctfsplit] <= j) {
					argdict[ctfsplit] = j + 1;
					printf("New CTF at %d\n", i);
				}
				defocus_val[j] = ctf->get_defocus();
				bfactor_val[j] = ctf->get_bfactor();
				outfile[strlen(outfile) - 4] = 0;
				sprintf(outfile + strlen(outfile), ".%02d.img", j);

				if( ld )
				{
					delete ld;
					ld = 0;
				}
				ld = d->copy();
			}

			float sigma = d->get_attr("sigma");
			if (!Util::goodf(&sigma)) {
				LOGWARN("Warning! bad Sigma for image %d in file '%s'", i, argv[1]);
				continue;
			}
#if 0
			if (argdict[norefs] && d->get_nimg() < 0) {
				continue;
			}
#endif
			if (argdict[edgenorm]) {
				d->process_inplace("normalize.circlemean");
			}

			if (norm == 1) {
				d->process_inplace("normalize");
			}
			else if (norm == 2) {
				(*d) *= nsig / sigma;
				(*d) += nmean - (float)d->get_attr("mean");
			}

			if (argdict[flip]) {
				d->process_inplace("xform.flip", Dict("axis", "y"));
			}

			if (argdict[invert]) {
				(*d) *= -1;
			}

			if (ramp) {
				d->process_inplace("eman1.filter.ramp", Dict("intercept", 1, "slope", ramp));	    
			}

			int y = d->get_ysize();

			if (argdict[setsfpairs]) {
				if (mask) {
					d->process_inplace(argfilters[noisemask], Dict("outer_radius", mask));
				}
				EMData *dataf = d->do_fft();

				float x0 = 0;
				float step = 0.5;
	    
				if (i % 2 == 0) {
					sfcurve1 = dataf->calc_radial_dist(nx, x0, step);
				}
				else {
					vector < float >sfcurve2 = dataf->calc_radial_dist(nx, x0, step);

					for (int j = 0; j < nx; j++) {
						if (sfcurve1[j] > 0 && sfcurve2[j] > 0) {
							sfcurve2[j] = sqrt(sfcurve1[j] / sfcurve2[j]);
						}
						else {
							sfcurve2[j] = 0;
						}
					}

					dataf->apply_radial_func(x0, step, sfcurve2);
					if( d )
					{
						delete d;
						d = 0;
					}
					d = dataf->do_ift();
				}
				
				if( dataf )
				{
					delete dataf;
					dataf = 0;
				}
			}

			float Xlp = lp;
			float Xtlp = tlp;
			float Xhp = hp;
			float Xsharphp = sharphp;

			if (apix > 0 && lp) {
				Xlp = y * apix / lp;
			}

			if (apix > 0 && tlp) {
				Xtlp = y * apix / tlp;
			}

			if (apix > 0 && hp) {
				Xhp = y * apix / hp;
			}

			if (apix > 0 && sharphp) {
				Xsharphp = y * apix / sharphp;
			}

			if (Xlp || Xhp) {
				d->process_inplace("filter.lowpass.gauss", Dict("cutoff_abs", Xhp == 0 ? -10.0 : Xhp));
				d->process_inplace("filter.highpass.tanh", Dict("cutoff_abs", Xlp == 0 ? 100000.0 : Xlp));
			}

			if (Xtlp) {
				d->process_inplace("filter.lowpass.tanh", Dict("cutoff_abs", -10.0));
				d->process_inplace("filter.highpass.tanh", Dict("cutoff_abs", Xtlp));
			}

			if (Xsharphp) {
				d->process_inplace("filter.lowpass.tophat", Dict("cutoff_abs", Xsharphp));
				d->process_inplace("filter.highpass.tophat", Dict("cutoff_abs", 100000.0));
			}

			if (mask) {
				d->process_inplace(argfilters[noisemask], Dict("outer_radius", mask));
			}

			if (imask > 0) {
				d->process_inplace("mask.sharp", Dict("inner_radius", imask, "value", 0));
			}

			if (automask) {
				d->process_inplace("mask.auto2d", Dict("threshold", automask));
			}

			// uses correlation with 180 deg rot for centering
			if (argdict[acfcenter]) {
				Dict params;
				params["intonly"] = 1;
				params["maxshift"] = d->get_xsize() / 4;
				d->align("translate", 0, params);
				//d->rotate_translate();
			}

			if (argdict[center]) {
				d->process_inplace("xform.centerofmass", Dict("int_shift_only", 1));
			}

			if (argdict[phot]) {
				d->process_inplace("xform.phaseorigin");
			}

			if (anoise) {
				d->process_inplace("math.addnoise");
			}

			if (argdict[rfp]) {
				EMData *e = d->make_rotational_footprint();
				e->append_image("rfp.hed");
			}


			if (rot || dx || dy || rize) {
				if (!rize) {
					d->rotate_translate(rot * M_PI / 180.0f, 0, 0, dx, dy, 0);
				}
				else {
					if (rizef && rand() % 2) {
						d->process_inplace("xform.flip", Dict("axis", "y"));
					}

					if (rizeda > 0) {
						d->rotate(Util::get_frand(-rizeda / 2.0f, rizeda / 2.0f), 0, 0);
					}
		
					if (rizedx > 0) {
						d->translate(Util::get_gauss_rand(0, rizedx),
									 Util::get_gauss_rand(0, rizedx), 0.0f);
					}
				}
			}

			Dict rr = d->get_transform().get_rotation("eman");
			//int nimg = d->get_nimg();

			if (scale < 1.0) {
				Transform t;
				t.set_scale(scale);
				d->transform(t);
			}

			if (clipx > 0) {
				EMData *e = d->get_clip(Region((d->get_xsize() - clipx) / 2,
											   (d->get_ysize() - clipy) / 2,
											   clipx, clipy));
				if( d )
				{
					delete d;
					d = 0;
				}
				d = e;
			}

			if (scale > 1.0) {
				Transform t;
				t.set_scale(scale);
				d->transform(t);
			}

			d->set_rotation((float)rr["alt"], (float)rr["az"], (float)rr["phi"]);
			//d->setNImg(nimg);

			if (fabs((float)shrink) > 1) {
				d->median_shrink((int) fabs((float)shrink));
			}

			if (argdict[rotav]) {
				d->process_inplace("math.rotationalaverage");
			}

			if (csym > 1) {
				EMData *f = d->copy();

				EMData *e = f->copy();
				for (int j = 1; j < csym; j++) {
					e->rotate(j * M_PI * 2.0f / csym, 0, 0);
					(*d) += (*e);
				}

				*d *= (1.0f / csym);

				if( e )
				{
					delete e;
					e = 0;
				}
				if( f )
				{
					delete f;
					f = 0;
				}
			}

			if (argdict[rsub]) {
				d->process_inplace("math.rotationalsubtract");
			}

			if (scl) {
				EMData *sc = new EMData();
				if (sclmd == 0) {
					sc->common_lines_real(d, d, scl, true);
				}
				else {
					EMData *e = d->copy();
					e->process_inplace("xform.phaseorigin");

					if (sclmd == 1) {
						sc->common_lines(e, e, sclmd, scl, true);
						sc->process_inplace("math.linear", Dict("shift", -90.0, "scale", -1.0));
					}
					else if (sclmd == 2) {
						sc->common_lines(e, e, sclmd, scl, true);
					}
					else {
						LOGERR("Invalid common-lines mode");
						exit(1);
					}
					if( e )
					{
						delete e;
						e = 0;
					}
				}
				if( d )
				{
					delete d;
					d = 0;
				}
					d = sc;
				
			}

			if (argdict[radon]) {
				EMData *r = d->do_radon();
				if( d )
				{
					delete d;
					d = 0;
				}
				d = r;
			}

			if (rfilt) {
				d->process_inplace(filtername, params_dict);
			}

			if (bliter > 0 && blwidth > 0) {
				Dict p;
				p["distance_sigma"] = blsigma1;
				p["value_sigma"] = blsigma2;
				p["niter"] = bliter;
				p["half_width"] = blwidth;

				d->process_inplace("bilateral", p);
			}

#if 0
			if (filefilt) {
				EMData *d2 = d->do_fft();
				if( d )
				{
					delete d;
					d = 0;
				}

				d2->apply_radial_func(nxyd, xd[0], xd[1] - xd[0], yd, true);
				d = d2->do_ift();
			}
#endif

			if (argdict[avgonly]) {
				if (!average) {
					average = d->copy();
				}
				else {
					(*average) += (*d);
				}
				continue;
			}
#if 0
			if (fftavgfile) {
				if (!fftavg) {
					fftavg = new EMData;
					fftavg->setSize(d->get_xsize() + 2, d->get_ysize());
					fftavg->setComplex(1);
					fftavg->zero();
				}
				d->applyMask(-1, 6);
				d->normalize();
				EMData *df = d->do_fft();
				df->multConst(df->get_ysize());	// just for scaling of the intensity level
				fftavg->addIncoherent(df);

				if( df )
				{
					delete df;
					df = 0;
				}
				continue;		// no writing yet
			}
#endif

#if 0
			if (strlen(sfout)) {
	    
				int y = d->get_ysize();
				float *curve = (float *) malloc(y * 4 * sizeof(float));
				EMData *dataf = d->do_fft();

				dataf->calc_radial_dist(y, 0, .5, curve);
				// single images go to outfile
				if (n1 == 0) {
					save_data(0, 1.0 / (apix * 2.0 * y), curve, y, sfout);
				}
				// multi-images get numbered
				else {
					sprintf(outfile2, "%s.%03d", sfout, i + 100);
					save_data(0, 1.0 / (apix * 2.0 * y), curve, y, outfile2);
				}
				if (sfout_n) {
					for (int j = 0; j < sfout_n; j++) {
						dataf->calc_radial_dist(y, 0, .5, curve, j * 2 * M_PI / sfout_n,
												2 * M_PI / sfout_n);
						// single images go to outfile
						if (n1 == 0) {
							sprintf(outfile2, "%s-%d-%d.pwr", file_basename(sfout), j, sfout_n);
							save_data(0, 1.0 / (apix * 2.0 * y), curve, y, outfile2);
						}
						// multi-images get numbered
						else {
							sprintf(outfile2, "%s-%d-%d.pwr.%03d", file_basename(sfout), j, sfout_n,
									i + 100);
							save_data(0, 1.0 / (apix * 2.0 * y), curve, y, outfile2);
						}
					}
				}
			}
#endif

			if (argdict[mrc]) {
				if (n1 != 0) {
					sprintf(outfile, "%03d.%s", i + 100, outfile);
				}
				d->write_image(outfile, 0, EMUtil::IMAGE_MRC);
			}
			else if (argdict[spidersingle]) {
				if (n1 != 0) {
					sprintf(spiderformat, "%s%%0%dd.spi", outfile, (int) (log10((float) i)) + 1);
					sprintf(outfile, spiderformat, i);
				}

				d->write_image(outfile, 0, EMUtil::IMAGE_SINGLE_SPIDER);
			}
			else if (argdict[hdf]) {
				d->append_image(outfile,  EMUtil::IMAGE_HDF);
			}
			else if (argdict[em]) {
				d->write_image(outfile, 0, EMUtil::IMAGE_EM);
			}
			else if (argdict[pif]) {
				if (n1 != 0) {
					sprintf(outfile, "%03d.%s", i + 100, outfile);
				}

				d->write_image(outfile, 0, EMUtil::IMAGE_PIF);
			}
			else if (argdict[png]) {
				if (n1 != 0) {
					sprintf(outfile, "%03d.%s", i + 100, outfile);
				}

				d->write_image(outfile, 0, EMUtil::IMAGE_PNG);
			}
			else if (pgm) {
				if (pgmlo >= pgmhi) {
					pgmlo = (int) d->get_attr("minimum");
					pgmhi = (int) d->get_attr("maximum");
				}

				d->set_attr("min_gray", pgmlo);
				d->set_attr("max_gray", pgmhi);

				if (n1 != 0) {
					sprintf(outfile, "%03d.%s", i + 100, outfile);
				}
				d->write_image(outfile, 0, EMUtil::IMAGE_PGM);
			}
			else if (argdict[spider]) {
				d->write_image(outfile, 0, EMUtil::IMAGE_SPIDER);
			}
			else if (argdict[vtk]) {
				d->write_image(outfile, 0, EMUtil::IMAGE_VTK);
			}
			else {
				if (argdict[inplace]) {
					d->write_image(outfile, i);
				}
				else if (argdict[mraprep]) {
					char nf[1024];
					sprintf(nf, "%s%04d.lst", outfile, i);
					d->write_image(nf, 0, EMUtil::IMAGE_LST);
				}
				else {
					d->append_image(outfile);
				}
			}

			if (interleaved_file) {
				d->read_image(interleaved_file, i);
				d->write_image(outfile, -1, EMUtil::IMAGE_IMAGIC);
			}

		}
	}
	catch(E2Exception &e ) {
		printf("%s\n", e.what());
	}
	
    if (average) {
		//average->setNImg(n1-n0+1);
		average->process_inplace("normalize");
		average->write_image(outfile, -1);
    }
#if 0
    if (fftavg) {
		fftavg->multConst(1. / sqrt((float) n1 - n0 + 1));
		fftavg->writeMRC(fftavgfile);
    }
#endif
    int n_outimg = EMUtil::get_image_count(argv[2]);
    printf("%d images\n", n_outimg);

    if (d) {
		delete d;
		d = 0;
    }

    if (ld) {
		delete ld;
		ld = 0;
    }

	if( defocus_val )
	{
    	delete [] defocus_val;
    	defocus_val = 0;
	}

	if( bfactor_val )
	{
    	delete [] bfactor_val;
    	bfactor_val = 0;
	}
    
    return 0;
}
