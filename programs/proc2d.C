#include "emdata.h"
#include "log.h"
#include <map>

using namespace EMAN;
using std::map;

void usage(const char *progname)
{
    printf("\n%s inputfile outputfile [-vN]\n", progname);
}

int main(int argc, char *argv[])
{
#if 0
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
    string plot = "plot";
    string avgonly = "average";
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

    string em = "em";
    string spider-single = "spider-single";
    string spiderswap-single = "spiderswap-single";
    string spiderswap = "spiderswap";

    string spider = "spider";
    string rotav = "rotav";
    string inplace = "inplace";
    string rsub = "rsub";

    string flip = "flip";
    string edgenorm = "edgenorm";
    string radon = "radon";
    
    
    map<string, bool> argdict;
    map<string, string> argfilters;
    
    float apix = 0;
    float lp = 0;
    float tlp = 0;

    int norm = 0;
    float nmean = 0;
    float nsig = 0;

    float ramp = 0;
    int mask = 0;
    
    int n0 = 0;
    int n1 = 0;

    argfilters[noisemask] = "MakeRadiusSquared";
    
    for (int i = 3; i < argc; i++) {
	const char *arg = argv[i];
	
	if (Util::get_str_float(arg, "apix=", &pix)) {
	}
	else if (Util::get_str_float(arg, "lp=", &lp)) {
	}
	else if (Util::get_str_float(arg, "tlp=", &tlp)) {
	}
	else if (Util::get_str_int(arg, "first=", &n0)) {
	}
	else if (Util::get_str_int(arg, "last=", &n1)) {
	}
	else if (Util::get_str_float(arg, "ramp=", &ramp)) {
	}
	else if (Util::get_str_int(arg, "mask=", &mask)) {
	}
	else if (Util::sstrncmp(arg, "norm")) {
	    if (arg[4] == '=') {
		norm = 2;
		sscanf(arg+5,"%f,%f", &nmean, &nsig);
	    }
	    else {
		norm = 1;
	    }
	}
	else if (Util::sstrncmp(arg, noisemask)) {
	    argfilters[noisemask] = "NoiseMask";
	}
	else {
	    argdict[argv[i]] = true;
	}
    }

    
    
    EMData *d = new EMData();
    
    int nimg = = EMUtil::get_image_count(imagefile);
    if (nimg <= n1 || n1 < 0) {
	n1 = nimg - 1;
    }
    const char *outfile = argv[2];
    EMData *ld = new EMData();
    
    for (int i = n0; i < n1; i++) {
	d->read_image(argv[1], i);
	
	if (ctfsplit && (i == n0 || !EMUtil::is_same_ctf(d, ld))) {
	    Ctf *ctf = d->get_ctf();

	    int j = 0;
	    for (j = 1; j < ctfsplit; j++) {
		if (defocus_val[j] == ctf->get_defocus() &&
		    bfactor_val[j] == ctf->get_bfactor()) {
		    break;
		}
	    }
	    if (ctfsplit <= j) {
		ctfsplit=j+1;
		printf("New CTF at %d\n",i);
	    }
	    dfval[j]=ctf->get_defocus();
	    bval[j]=ctf->get_bfactor();
	    outfile[strlen(outfile)-4]=0;
	    sprintf(outfile+strlen(outfile),".%02d.img",j);
	    
	    delete ld;
	    ld=d->copy(false,false);
	}

	float sigma = d->get_sigma();
	if (goodf(sigma)) {
	    continue;
	}
	if (sigma == 0) {
	    Log::logger()->warn("Warning! Sigma=0 for image %d", i);
	}
	
	if (argdict[norefs] && d->get_nimg<0) {
	    continue;
	}
	
	if (argdict[edgenorm]) {
	    d->filter("CircleMeanNormalize");
	}
	
	if (norm==1) {
	    d->filter("StdNormalize");
	}
	else if (norm==2) {
	    (*d) *= nsig/sigma;
	    (*d) += nmean - d->get_mean();
	}
	
	if (argdict[flip]) {
	    d->vertical_flip();
	}
	
	if (argdict[invert]) {
	    (*d) *= -1;
	}

	if (ramp) {
	    d->filter("LinearRamp", Dict("intercept", EMObject(1), "slope", EMObject(ramp)));
	    d->filter(ramp,1.,6);
	}
	int y=d->get_ysize();

	if (argdict[setsfpairs]) {
	    if (mask) {
		d->filter(argfilters[noisemask], Dict("outer_radius", EMObject(mask)));
	    }
	    EMData *dataf=d->doFFT();
	    d->gimmeFFT();
	    if (i%2==0) {
		dataf->calcRadDist(nx,0,.5,sfcurve1);
	    }
	    else {
		dataf->calcRadDist(nx,0,.5,sfcurve2);
		for (int j=0; j<nx; j++) {
		    if (sfcurve1[j]>0 && sfcurve2[j]>0) sfcurve2[j]=sqrt(sfcurve1[j]/sfcurve2[j]);
		    else sfcurve2[j]=0;
		}

		dataf->applyRadFn(nx,0,.5,sfcurve2);
		delete d;
		d=dataf->doIFT();
		dataf->gimmeFFT();
	    }
	    delete dataf;
	}

	if (snrfilt) {
	    float *ctf;
	    EMData *d2;
	    if (wiener) ctf=d->ctfCurve(6,&SF);
	    else ctf=d->ctfCurve(4,&SF);
	    d->edgeNormalize();
	    EMData *d3=d->clip(-d->xSize()/2,-d->ySize()/2,d->xSize()*2,d->ySize()*2);
	    d2=d3->doFFT();
	    d3->gimmeFFT();
	    d2->applyRadFn(CTFOS*d->ySize()/2,0,2.0/CTFOS,ctf,0);
	    delete d3;
	    delete d;
	    d3=d2->doIFT();
	    d=d3->clip(d3->xSize()/4,d3->ySize()/4,d3->xSize()/2,d3->ySize()/2);
	    d->setParent(NULL);
	    delete d3;
	    d2->gimmeFFT();
	    delete d2;
	    free(ctf);
	}

	float Xlp,Xhp,Xhp2,Xtlp;
	if (apix>0 && lp) Xlp=y*apix/lp; else Xlp=lp;
	if (apix>0 && tlp) Xtlp=y*apix/tlp; else Xtlp=tlp;
	if (apix>0 && hp) Xhp=y*apix/hp; else Xhp=hp;
	if (apix>0 && hp2) Xhp2=y*apix/hp2; else Xhp2=hp2;

	if (Xlp || Xhp) {  d->filter(Xhp==0?-10.0:Xhp,Xlp==0?100000.0:Xlp,1); }
	if (Xtlp) d->filter(-10.0,Xtlp,2);
	if (Xhp2) d->filter(Xhp2,100000.0,0);

	if (mask) d->applyMask(mask,noisemask);
	if (imask>0) d->applyMask(imask,5);

	if (amask) d->autoMask(amask,.1);

	// uses correlation with 180 deg rot for centering
	if (acfcenter) { d->transAlign(NULL,0,1,d->xSize()/4); d->rotateAndTranslate(); }

	// uses center of mass to center
	if (center) d->cmCenter(1);
	if (phot) d->toCorner();
	//	if (d->Sigma()==0 || !finite(d->Sigma())) continue;

	if (anoise) {
	    dat=d->getData();
	    for (j=0; j<d->xSize()*d->ySize(); j++) dat[j]+=grand(anoise,anoise/2);
	    d->doneData();
	}

	if (rfp) { 
	    e=d->makeRFP(); 
	    e->writeImage("rfp.hed",-1); 
	}

	if (rot||dx||dy||rize) {
	    Euler oldang=Euler(*d->getEuler());
	    float dx0=d->Dx(),	dy0=d->Dy();
	    d->setParent(NULL);
	    d->setRAlign(rot*PI/180.0,0,0);
	    d->setTAlign(dx,dy,0);
	    if (rize) {
		if (rizeda>0) d->setRAlign(frand(-rizeda/2.0,rizeda/2.0));
		if (rizedx>0) d->setTAlign(grand(0,rizedx),grand(0,rizedx),0);
		if (rizef && random()%2) d->vFlip();
	    }
	    d->rotateAndTranslate();
	    d->setRAlign(oldang.alt(),oldang.az(),oldang.phi()+rot*PI/180.0);
	    d->setTAlign(dx0+dx,dy0+dy,0);
	}

	Euler ee = *d->getEuler();
	int ni = d->NImg();
	if (scale<1.0) { d->setRAlign(0,0,0); d->setTAlign(0,0,0); d->rotateAndTranslate(scale); }
	if (clip>0) {
	    e=d->clip((d->xSize()-clip)/2,(d->ySize()-clipy)/2,clip,clipy);
	    delete d;
	    d=e;
	}
	if (scale>1.0) { d->setRAlign(0,0,0); d->setTAlign(0,0,0); d->rotateAndTranslate(scale); }
	d->setRAlign(&ee);
	d->setNImg(ni);

	if (shrink>1) d->medianShrink(shrink);
	else if (shrink<-1) d->meanShrink(-shrink);

	// this is handled by the library now
	/*	if (abs(shrink)>1 || scale!=1.0) {
		if (d->hasCTF()) {
		float *ctf=d->getCTF();
		ctf[10]*=fabs((float)abs(shrink))/scale;
		}
		}*/

	if (rotav) {
	    d->radialAverage();
	}


	if (csym>1) {
	    f=d->copy();
	    f->setParent(NULL);
	    e=f->copy();
	    for (j=1; j<csym; j++) {
		e->setRAlign(j*PI*2.0/(float)csym,0,0);
		e->rotateAndTranslate();
		d->add(e);
	    }
	    d->multConst(1.0/(float)csym);
	    delete e;
	    delete f;
	}

	if (rsub) {
	    d->radialSubtract();
	}

	if (scl) {
	    EMData *sc = new EMData;
	    if (sclmd==0) sc->commonLinesR(d,d,scl,1);
	    else {
		//			EMData *e=d->clip(-d->xSize()/2,-d->ySize()/2,d->xSize()*2,d->ySize()*2);
		EMData *e=d->copy();
		e->toCorner();
		if (sclmd==1) {
		    sc->commonLines(e,e,1,scl,1);
		    sc->realFilter(9,-90.0,-1.0);
		}
		else if (sclmd==2) {
		    sc->commonLines(e,e,2,scl,1);
		}
		else error(ERR_EXIT,"Invalid common-lines mode");
		delete e;
	    }
	    delete d;
	    d=sc;
	}

	if(radon) {
	    EMData *r=d->doRadon();
	    delete d;
	    d=r;
	}


	if (rfilt>=0) d->realFilter(rfilt,rft1,rft2,rft3);

	if (bliter >0 && blwidth>0) {	// bilateral filtering
	    d->BilateralFilter(blsigma1, blsigma2, bliter, blwidth);
	}
	// This filters an image using a filter specified as a 2 column text file with s and the filter amp
	if (filefilt) {
	    EMData *d2;
	    d2=d->doFFT();
	    d->gimmeFFT();
	    delete d;
	    d2->applyRadFn(nxyd,xd[0],xd[1]-xd[0],yd,1);
	    d=d2->doIFT();
	    d2->gimmeFFT();
	}

	if (avgonly) {
	    if (!average) average=d->copy(0,0);
	    else average->add(d);
	    continue;
	}

	if (fftavgfile) {
	    if (!fftavg) { fftavg=new EMData; fftavg->setSize(d->xSize()+2,d->ySize());	fftavg->setComplex(1);	fftavg->zero(); }
	    d->applyMask(-1,6);
	    d->normalize();
	    EMData* df=d->doFFT();
	    df->multConst(df->ySize());	// just for scaling of the intensity level
	    fftavg->addIncoherent(df);	
	    d->gimmeFFT();
	    delete df;
	    continue;	// no writing yet
	}
	

	if (strlen(sfout)) {
	    char outfile2[1024];
	    int y = d->ySize();
	    float *curve = (float *)malloc(y*4*sizeof(float));
	    EMData *dataf=d->doFFT();
	    d->gimmeFFT();
	    dataf->calcRadDist(y,0,.5,curve);
	    // single images go to outfile
	    if (n1==0) {
		save_data(0,1.0/(apix*2.0*y),curve,y,sfout);
	    }
	    // multi-images get numbered
	    else {
		sprintf(outfile2,"%s.%03d",sfout,i+100);
		save_data(0,1.0/(apix*2.0*y),curve,y,outfile2);
	    }
	    if (sfout_n) {
		for(int j=0; j<sfout_n; j++){
		    dataf->calcRadDist(y,0,.5,curve,j*2*PI/sfout_n,2*PI/sfout_n);
		    // single images go to outfile
		    if (n1==0) {
			sprintf(outfile2,"%s-%d-%d.pwr",file_basename(sfout),j,sfout_n);
			save_data(0,1.0/(apix*2.0*y),curve,y,outfile2);
		    }
		    // multi-images get numbered
		    else {
			sprintf(outfile2,"%s-%d-%d.pwr.%03d",file_basename(sfout),j,sfout_n,i+100);
			save_data(0,1.0/(apix*2.0*y),curve,y,outfile2);
		    }
		}
	    }
	}

	if (mrc) {
	    // single images go to outfile
	    if (n1==0) {
		d->writeMRC(outfile,mrc==2?0:2);
	    }
	    // multi-images get numbered
	    else {
		sprintf(outfile,"%03d.%s",i+100,argv[2]);
		d->writeMRC(outfile,mrc==2?0:2);
	    }
	}
	else if (spidersingle) {
	    if (n1==0) {
		d->writeSingleSPIDER(outfile,spidersingle-1);
	    }
	    else {
		sprintf(spiderfile,spiderformat,i);
		d->writeSingleSPIDER(spiderfile,spidersingle-1);
	    }
	}
	else if (hdf) {
	    d->writeHDF(outfile, -1);
	}
	else if (em) {
	    d->writeEM(outfile);
	}
	else if (pif) {
	    if (n1==0) {
		d->writePIF(outfile, false);
	    }
	    else {
		sprintf(outfile,"%03d.%s",i+100,argv[2]);
		d->writePIF(outfile, false);
	    }
	}
	else if (png) {
	    if (n1==0) {
		d->writePNG(outfile);
	    }
	    else {
		sprintf(outfile,"%03d.%s",i+100,argv[2]);
		d->writePNG(outfile);
	    }
	}
	else if (pgm) {
	    // single images go to outfile
	    if (n1==0) {
		if (pgmlo>=pgmhi) d->writePGM(outfile,d->Min(),d->Max());
		else d->writePGM(outfile,pgmlo,pgmhi);
	    }
	    // multi-images get numbered
	    else {
		sprintf(outfile,"%03d.%s",i+100,argv[2]);
		if (pgmlo>=pgmhi) d->writePGM(outfile,d->Min(),d->Max());
		else d->writePGM(outfile,pgmlo,pgmhi);
	    }
	}
	else if (spider) {
	    if (spider==1) d->writeSPIDER(outfile,-1);
	    else d->writeSPIDER(outfile,-1,1);	// write byte-swapped
	}
	else {
	    if (inpl) d->writeImage(outfile,i);
	    else if (mraprep) {
		char nf[1024];
		sprintf(nf,"%s%04d.lst",outfile,i);
		d->writeLST(nf,0);
	    }
	    else d->writeImage(outfile,-1);
	}

	// interleaved images not processed
	if (il) {
	    d->readImage(il,i);
	    d->writeImage(outfile,-1);
	}

	//	if (rfp) delete e;	// auotmatically freed by emdata
    }

    if (average) {
	average->setNImg(n1-n0+1);
	average->normalize();
	average->writeImage(outfile,-1);
    }

    if (fftavg) {
	fftavg->multConst(1./sqrt((float)n1-n0+1));
	fftavg->writeMRC(fftavgfile);
    }

    LOG(Ref,argv[2],LOG_OUTFILE,NULL);
    LOG(Ref,NULL,LOG_END,NULL);

    printf("%d images\n",fileCount(argv[2]));
}


    
#endif

return 0;
}
