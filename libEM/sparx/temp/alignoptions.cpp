#include "alignoptions.h"

AlignOptions::AlignOptions()
{
    mask3D     = NULL;
    first_ring = 1;
    last_ring  = 26;
    rstep      = 1;
    ri         = 26;
    xrng       = 1.0;
    yrng       = 1.0;
    step       = 1.0;
    dtheta     = 15.0;
    snr        = 1.0;
    symmetry   = "c1";
    CTF        = false;
}

AlignOptions::AlignOptions(const AlignOptions& old_options)
{
    mask3D     = old_options.mask3D;
    first_ring = old_options.first_ring;
    last_ring  = old_options.last_ring;
    rstep      = old_options.rstep;
    ri         = old_options.ri;
    xrng       = old_options.xrng;
    yrng       = old_options.yrng;
    step       = old_options.step;
    dtheta     = old_options.dtheta;
    snr        = old_options.snr;
    symmetry   = old_options.symmetry;
    CTF        = old_options.CTF;
}

AlignOptions::~AlignOptions()
{
}

void AlignOptions::set_mask3D(EMData * in_mask3D){
    mask3D = in_mask3D;
}
void AlignOptions::set_first_ring(int in_first_ring){
    first_ring = in_first_ring;
}
void AlignOptions::set_last_ring(int in_last_ring){
    last_ring = in_last_ring;
}
void AlignOptions::set_rstep(int in_rstep){
    rstep = in_rstep;
}
void AlignOptions::set_ri(int in_ri){
    ri = in_ri;
}
void AlignOptions::set_xrng(float in_xrng){
    xrng = in_xrng;
}
void AlignOptions::set_yrng(float in_yrng){
    yrng = in_yrng;
}
void AlignOptions::set_step(float in_step){
    step = in_step;
}
void AlignOptions::set_dtheta(float in_dtheta){
    dtheta = in_dtheta;
}
void AlignOptions::set_snr(float in_snr){
    snr = in_snr;
}
void AlignOptions::set_symmetry(std::string in_symmetry){
    symmetry = in_symmetry;
}
void AlignOptions::set_CTF(bool in_CTF){
    CTF = in_CTF;
}

EMData * AlignOptions::get_mask3D(){
    return mask3D;
}
int   AlignOptions::get_first_ring(){
    return first_ring;
}
int   AlignOptions::get_last_ring(){
    return last_ring;
}
int   AlignOptions::get_rstep(){
    return rstep;
}
int   AlignOptions::get_ri(){
    return ri;
}
float AlignOptions::get_xrng(){
    return xrng;
}
float AlignOptions::get_yrng(){
    return yrng;
}
float AlignOptions::get_step(){
    return step;
}
float AlignOptions::get_dtheta(){
    return dtheta;
}
float AlignOptions::get_snr(){
    return snr;
}
std::string AlignOptions::get_symmetry(){
    return symmetry;
}
bool AlignOptions::get_CTF(){
    return CTF;
}
