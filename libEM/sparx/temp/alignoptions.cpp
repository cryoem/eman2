#include "alignoptions.h"

AlignOptions::AlignOptions()
{
    first_ring = 1;
    last_ring  = 26;
    rstep      = 1;
    ri         = 26;
    xrng       = 1.0;
    yrng       = 1.0;
    step       = 1.0;
    dtheta     = 15.0;
}

AlignOptions::AlignOptions(const AlignOptions& old_options)
{
    first_ring = old_options.first_ring;
    last_ring  = old_options.last_ring;
    rstep      = old_options.rstep;
    ri         = old_options.ri;
    xrng       = old_options.xrng;
    yrng       = old_options.yrng;
    step       = old_options.step;
    dtheta     = old_options.dtheta;
    
}

AlignOptions::~AlignOptions()
{
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
