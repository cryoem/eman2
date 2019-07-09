/*
 * Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
 * Please do not copy or modify this file without written consent of the author.
 * Copyright (c) 2000-2019 The University of Texas - Houston Medical School
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
 */
#ifndef ALIGNOPTIONS_H
#define ALIGNOPTIONS_H

#include "mpi.h"
#include "emdata.h"
#include <string>

using namespace EMAN;

class AlignOptions {
public:
    AlignOptions();
    AlignOptions(Vec3i& volsize);
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
    void set_have_angles(bool have_angles);
    void set_ref_angle_type(std::string ref_angle_type);
    void set_use_sirt(bool use_sirt);
    void set_sirt_tol(float sirt_tol);
    void set_sirt_lam(float sirt_lam);
    void set_sirt_maxit(int sirt_maxit);
    void set_maxit(int maxit);

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
    bool  get_have_angles();
    std::string get_ref_angle_type();
    bool  get_use_sirt();
    float get_sirt_tol();
    float get_sirt_lam();
    int   get_sirt_maxit();
    int   get_maxit();
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
    bool have_angles;
    std::string ref_angle_type;
    bool use_sirt;
    float sirt_tol;
    float sirt_lam;
    int sirt_maxit;
    int maxit;
};
#endif // ALIGNOPTIONS_H
