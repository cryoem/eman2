#ifndef ALIGNOPTIONS_H
#define ALIGNOPTIONS_H

#include "emdata.h"
#include <string>

using namespace EMAN;

class AlignOptions {
public:
    AlignOptions();
    AlignOptions(const AlignOptions& old_options);
    virtual ~AlignOptions();

    void set_mask3D(EMData * mask3D);
    void set_first_ring(int first_ring);
    void set_last_ring(int last_ring);
    void set_rstep(int rstep);
    void set_ri(int ri);
    void set_xrng(float xrng);
    void set_yrng(float yrng);
    void set_step(float step);
    void set_dtheta(float dtheta);
    void set_snr(float snr);
    void set_symmetry(std::string symmetry);
    void set_CTF(bool CTF);

    EMData * get_mask3D();
    int   get_first_ring();
    int   get_last_ring();
    int   get_rstep();
    int   get_ri();
    float get_xrng();
    float get_yrng();
    float get_step();
    float get_dtheta();
    float get_snr();
    std::string get_symmetry();
    bool  get_CTF();

private:

    EMData * mask3D;
    int first_ring;
    int last_ring;
    int rstep;
    int ri;
    float xrng;
    float yrng;
    float step;
    float dtheta;
    float snr;
    std::string symmetry;
    bool CTF;

};
#endif // ALIGNOPTIONS_H
