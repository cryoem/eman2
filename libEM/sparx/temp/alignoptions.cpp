#include "alignoptions.h"

AlignOptions::AlignOptions()
{
    mask3D     = NULL;
    first_ring = 1;
    last_ring  = 0;
    rstep      = 1;
    ri         = 0;
    xrng       = 1.0;
    yrng       = 1.0;
    step       = 1.0;
    dtheta     = 15.0;
    snr        = 1.0;
    symmetry   = "c1";
    CTF        = false;
    have_angles = false;
    ref_angle_type = "P";
    use_sirt = true;
    sirt_tol = 1.0e-1;
    sirt_lam = 1.0e-4;
    sirt_maxit = 100;
    maxit      = 1;
}

AlignOptions::AlignOptions(Vec3i& volsize)
{
    int min_dim = (volsize[0] < volsize[1] ? volsize[0] : volsize[1]);
    min_dim = (min_dim < volsize[2] ? min_dim : volsize[2]);
    // min_dim = min(nx,ny,nz)
    mask3D     = NULL;
    first_ring = 1;
    last_ring  = min_dim / 2 - 2;
    rstep      = 1;
    ri         = min_dim / 2 - 2;
    xrng       = 1.0;
    yrng       = 1.0;
    step       = 1.0;
    dtheta     = 15.0;
    snr        = 1.0;
    symmetry   = "c1";
    CTF        = false;
    have_angles = false;
    ref_angle_type = "P";
    use_sirt = true;
    sirt_tol = 1.0e-1;
    sirt_lam = 1.0e-4;
    sirt_maxit = 100;
    maxit      = 1;
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
    have_angles = old_options.have_angles;
    ref_angle_type = old_options.ref_angle_type;
    use_sirt = old_options.use_sirt;
    sirt_tol = old_options.sirt_tol;
    sirt_lam = old_options.sirt_lam;
    sirt_maxit = old_options.sirt_maxit;
    maxit      = old_options.maxit;
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
void AlignOptions::set_have_angles(bool in_have_angles){
    have_angles = in_have_angles;
}
void AlignOptions::set_ref_angle_type(std::string in_ref_angle_type){
    ref_angle_type = in_ref_angle_type;
}
void AlignOptions::set_use_sirt(bool in_use_sirt){
    use_sirt = in_use_sirt;
}
void AlignOptions::set_sirt_tol(float in_sirt_tol){
    sirt_tol = in_sirt_tol;
}
void AlignOptions::set_sirt_lam(float in_sirt_lam){
    sirt_lam = in_sirt_lam;
}
void AlignOptions::set_sirt_maxit(int in_sirt_maxit){
    sirt_maxit = in_sirt_maxit;
}
void AlignOptions::set_maxit(int in_maxit){
    maxit = in_maxit;
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
bool AlignOptions::get_have_angles(){
    return have_angles;
}
std::string AlignOptions::get_ref_angle_type(){
    return ref_angle_type;
}
bool AlignOptions::get_use_sirt(){
    return use_sirt;
}
float AlignOptions::get_sirt_tol(){
    return sirt_tol;
}
float AlignOptions::get_sirt_lam(){
    return sirt_lam;
}
int AlignOptions::get_sirt_maxit(){
    return sirt_maxit;
}
int AlignOptions::get_maxit(){
    return maxit;
}
