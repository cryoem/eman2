#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from past.utils import old_div
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
# Author: Pawel Penczek, 4/4/2007 (Pawel.A.Penczek@uth.tmc.edu)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#

#  Implemented by PAP 03/02/09.  Based on Agarwal Acta Cryst. 1978, A34, 791.

# PDB sample line
#           1         2         3         4         5         6         7
# 01234567890123456789012345678901234567890123456789012345678901234567890123456789
# ATOM      1  N   ASP L   1      57.169  62.936   7.756  1.00 30.65      1ACY 129
# HEADER    COMPLEX(ANTIBODY/HIV-1 FRAGMENT)        10-FEB-94   1ACY      1ACY   2


import EMAN2_cppwrap
import EMAN2_meta
import math
import optparse
from ..libpy import sp_global_def
from ..libpy import sp_utilities
import sys
from builtins import range

# "2 DIGIT ATOM IN CAPSLOG": (atomic number, atomic weight),
atomdefs = {
    "H": (1.0, 1.00794),
    "C": (6.0, 12.0107),
    "A": (7.0, 14.00674),
    "N": (7.0, 14.00674),
    "O": (8.0, 15.9994),
    "MG": (12.0, 24.305),
    "P": (15.0, 30.973761),
    "S": (16.0, 32.066),
    "W": (18.0, 1.00794 * 2.0 + 15.9994),
    "AU": (79.0, 196.96655),
}


def run():
    progname = optparse.os.path.basename(sys.argv[0])
    usage = """%prog [options] input.pdb output.hdf

Converts a pdb file into an electron density map. 0,0,0 in PDB space will 
map to the center of the volume."""

    parser = optparse.OptionParser(usage=usage, version=EMAN2_meta.EMANVERSION)

    parser.add_option("--apix", "-A", type="float", help="Angstrom/voxel", default=1.0)
    parser.add_option(
        "--box", "-B", type="string", help="Box size in pixels, <xyz> or <x,y,z>"
    )
    parser.add_option(
        "--het", action="store_true", help="Include HET atoms in the map", default=False
    )
    parser.add_option(
        "--chains",
        type="string",
        help="String list of chain identifiers to include, e.g. 'ABEFG'; default: include all chains",
        default="",
    )
    parser.add_option(
        "--center",
        type="string",
        default="a",
        help="center: c - coordinates; a (default) - center of gravity; <x,y,z> - vector (in Angstrom) to subtract from all coordinates; n - none",
    )
    parser.add_option(
        "--O", action="store_true", default=False, help="use O system of coordinates"
    )
    parser.add_option(
        "--quiet", action="store_true", default=False, help="Verbose is the default"
    )
    parser.add_option(
        "--tr0",
        type="string",
        default="none",
        help="Filename of initial 3x4 transformation matrix",
    )
    parser.add_option(
        "--set_apix_value",
        action="store_true",
        help="Set apix value in header of the ouput map",
        default=False,
    )

    (options, args) = parser.parse_args()  #

    if len(args) < 2:
        sp_global_def.ERROR("Input and output files required")
        return

    if sp_global_def.CACHE_DISABLE:
        sp_utilities.disable_bdb_cache()

    chains = options.chains
    if chains == "":
        chains = None

    try:
        infile = open(args[0], "r")
    except:
        sp_global_def.ERROR("Cannot open input file")
        return

    aavg = [0, 0, 0]  # calculate atomic center
    asig = [0, 0, 0]  # to compute radius of gyration
    natm = 0
    atoms = []  # we store a list of atoms to process to avoid multiple passes
    nelec = 0
    mass = 0

    # read in initial-transformation file:
    if options.tr0 != "none":
        cols = sp_utilities.read_text_file(options.tr0, -1)
        txlist = []
        for i in range(3):
            txlist.append(cols[0][i])
            txlist.append(cols[1][i])
            txlist.append(cols[2][i])
            txlist.append(cols[3][i])
        tr0 = EMAN2_cppwrap.Transform(txlist)

    # parse the pdb file and pull out relevant atoms
    for line in infile:
        if line[:4] == "ATOM" or (line[:6] == "HETATM" and options.het):
            if chains and (line[21] not in chains):
                continue
            try:
                a = line[12:14].strip()
                aseq = int(line[6:11].strip())
                res = int(line[22:26].strip())
                if options.O:
                    x = float(line[38:46])
                    y = float(line[30:38])
                    z = -float(line[46:54])
                else:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
            except:
                sp_global_def.sxprint(
                    "PDB Parse error:\n%s\n'%s','%s','%s'  '%s','%s','%s'\n"
                    % (
                        line,
                        line[12:14],
                        line[6:11],
                        line[22:26],
                        line[30:38],
                        line[38:46],
                        line[46:54],
                    )
                )
                sp_global_def.sxprint(a, aseq, res, x, y, z)

            try:
                nelec += atomdefs[a.upper()][0]
                mass += atomdefs[a.upper()][1]
            except:
                sp_global_def.sxprint(("Unknown atom %s ignored at %d. You can modify the atomdefs info at the top of this file with the atomic number and atomic weight." % (a, aseq)))
                continue

            atoms.append([a, x, y, z])
            natm += 1

            if options.center == "a":
                aavg[0] += x * atomdefs[a.upper()][1]
                aavg[1] += y * atomdefs[a.upper()][1]
                aavg[2] += z * atomdefs[a.upper()][1]
                asig[0] += x ** 2 * atomdefs[a.upper()][1]
                asig[1] += y ** 2 * atomdefs[a.upper()][1]
                asig[2] += z ** 2 * atomdefs[a.upper()][1]
            else:
                aavg[0] += x
                aavg[1] += y
                aavg[2] += z
                asig[0] += x ** 2
                asig[1] += y ** 2
                asig[2] += z ** 2
    infile.close()

    if options.center == "a":
        rad_gyr = math.sqrt(
            old_div((asig[0] + asig[1] + asig[2]), mass)
            - (old_div(aavg[0], mass)) ** 2
            - (old_div(aavg[1], mass)) ** 2
            - (old_div(aavg[2], mass)) ** 2
        )
    else:
        rad_gyr = math.sqrt(
            old_div((asig[0] + asig[1] + asig[2]), natm)
            - (old_div(aavg[0], natm)) ** 2
            - (old_div(aavg[1], natm)) ** 2
            - (old_div(aavg[2], natm)) ** 2
        )

    if not options.quiet:
        sp_global_def.sxprint(
            "%d atoms; total charge = %d e-; mol mass = %.2f kDa; radius of gyration = %.2f A"
            % (natm, nelec, old_div(mass, 1000.0), rad_gyr)
        )

    # center PDB according to option:
    if options.center == "a":
        if not options.quiet:
            sp_global_def.sxprint(
                "center of gravity at %1.1f,%1.1f,%1.1f (center of volume at 0,0,0)"
                % (
                    old_div(aavg[0], mass),
                    old_div(aavg[1], mass),
                    old_div(aavg[2], mass),
                )
            )
        for i in range(len(atoms)):
            atoms[i][1] -= old_div(aavg[0], mass)
            atoms[i][2] -= old_div(aavg[1], mass)
            atoms[i][3] -= old_div(aavg[2], mass)
    if options.center == "c":
        if not options.quiet:
            sp_global_def.sxprint(
                "atomic center at %1.1f,%1.1f,%1.1f (center of volume at 0,0,0)"
                % (
                    old_div(aavg[0], natm),
                    old_div(aavg[1], natm),
                    old_div(aavg[2], natm),
                )
            )
        for i in range(len(atoms)):
            atoms[i][1] -= old_div(aavg[0], natm)
            atoms[i][2] -= old_div(aavg[1], natm)
            atoms[i][3] -= old_div(aavg[2], natm)
    spl = options.center.split(",")
    if len(spl) == 3:  # substract the given vector from all coordinates
        if not options.quiet:
            sp_global_def.sxprint(
                "vector to substract: %1.1f,%1.1f,%1.1f (center of volume at 0,0,0)"
                % (float(spl[0]), float(spl[1]), float(spl[2]))
            )
        for i in range(len(atoms)):
            atoms[i][1] -= float(spl[0])
            atoms[i][2] -= float(spl[1])
            atoms[i][3] -= float(spl[2])

    # apply given initial transformation (this used to be done before centering,
    # thereby loosing the translation. This is the right place to apply tr0):
    if options.tr0 != "none":
        if not options.quiet:
            sp_global_def.sxprint(
                "Applying initial transformation to PDB coordinates... "
            )
        for i in range(len(atoms)):
            atom_coords = EMAN2_cppwrap.Vec3f(atoms[i][1], atoms[i][2], atoms[i][3])
            new_atom_coords = tr0 * atom_coords
            atoms[i][1] = new_atom_coords[0]
            atoms[i][2] = new_atom_coords[1]
            atoms[i][3] = new_atom_coords[2]
        if not options.quiet:
            sp_global_def.sxprint("done.\n")

    # bounding box:
    amin = [atoms[0][1], atoms[0][2], atoms[0][3]]
    amax = [atoms[0][1], atoms[0][2], atoms[0][3]]
    for i in range(1, len(atoms)):
        for k in range(3):
            amin[k] = min(atoms[i][k + 1], amin[k])
            amax[k] = max(atoms[i][k + 1], amax[k])

    if not options.quiet:
        sp_global_def.sxprint(
            "Range of coordinates [A]: x: %7.2f - %7.2f" % (amin[0], amax[0])
        )
        sp_global_def.sxprint(
            "                          y: %7.2f - %7.2f" % (amin[1], amax[1])
        )
        sp_global_def.sxprint(
            "                          z: %7.2f - %7.2f" % (amin[2], amax[2])
        )

    # find the output box size, either user specified or from bounding box
    box = [0, 0, 0]
    try:
        spl = options.box.split(",")
        if len(spl) == 1:
            box[0] = box[1] = box[2] = int(spl[0])
        else:
            box[0] = int(spl[0])
            box[1] = int(spl[1])
            box[2] = int(spl[2])
    except:
        for i in range(3):
            box[i] = int(
                old_div(2 * max(math.fabs(amax[i]), math.fabs(amin[i])), options.apix)
            )
            #  Increase the box size by 1/4.
            box[i] += old_div(box[i], 4)

    if not options.quiet:
        sp_global_def.sxprint("Bounding box [pixels]: x: %5d " % box[0])
        sp_global_def.sxprint("                       y: %5d " % box[1])
        sp_global_def.sxprint("                       z: %5d " % box[2])

    # figure oversampled box size
    # bigb = max(box[0],box[1],box[2])
    fcbig = 1
    """Multiline Comment0"""
    if not options.quiet:
        sp_global_def.sxprint(
            "Box size: %d x %d x %d" % (box[0], box[1], box[2]),
            ",  oversampling ",
            fcbig,
        )

    # Calculate working dimensions
    pixelbig = old_div(options.apix, fcbig)
    bigbox = []
    for i in range(3):
        bigbox.append(box[i] * fcbig)

    # initialize the final output volume
    outmap = EMAN2_cppwrap.EMData(bigbox[0], bigbox[1], bigbox[2], True)
    nc = []
    for i in range(3):
        nc.append(old_div(bigbox[i], 2))
    # fill in the atoms
    for i in range(len(atoms)):
        # print "Adding %d '%s'"%(i,atoms[i][0])
        if not options.quiet and i % 1000 == 0:
            sp_global_def.sxprint("\r   %d" % i, end=" ")
            sys.stdout.flush()
        try:
            elec = atomdefs[atoms[i][0].upper()][0]
            # outmap[int(atoms[i][1]/pixelbig+bigbox[0]//2),int(atoms[i][2]/pixelbig+bigbox[1]//2),int(atoms[i][3]/pixelbig+bigbox[2]//2)] += elec
            for k in range(2):
                pz = old_div(atoms[i][3], pixelbig) + nc[2]
                dz = pz - int(pz)
                uz = ((1 - k) + (2 * k - 1) * dz) * elec
                for l in range(2):
                    py = old_div(atoms[i][2], pixelbig) + nc[1]
                    dy = py - int(py)
                    uy = ((1 - l) + (2 * l - 1) * dy) * uz
                    for m in range(2):
                        px = old_div(atoms[i][1], pixelbig) + nc[0]
                        dx = px - int(px)
                        outmap[int(px) + m, int(py) + l, int(pz) + k] += (
                            (1 - m) + (2 * m - 1) * dx
                        ) * uy
        except:
            sp_global_def.sxprint("Skipping %d '%s'" % (i, atoms[i][0]))

    if not options.quiet:
        sp_global_def.sxprint(
            "\r   %d\nConversion complete." % len(atoms)
        )  # ,"    Now shape atoms."
    """Multiline Comment1"""
    (filename_path, filextension) = optparse.os.path.splitext(args[1])
    if filextension == ".hdf":
        if options.set_apix_value:
            outmap.set_attr("apix_x", options.apix)
            outmap.set_attr("apix_y", options.apix)
            outmap.set_attr("apix_z", options.apix)
            outmap.set_attr("pixel_size", options.apix)
        else:
            sp_global_def.sxprint("Pixel_size is not set in the header!")

        outmap.write_image(args[1], 0, EMAN2_cppwrap.EMUtil.ImageType.IMAGE_HDF)

    elif filextension == ".spi":
        outmap.write_image(
            args[1], 0, EMAN2_cppwrap.EMUtil.ImageType.IMAGE_SINGLE_SPIDER
        )

    else:
        sp_global_def.ERROR("Unknown image type")
        return

def main():
    sp_global_def.print_timestamp("Start")
    sp_global_def.write_command()
    run()
    sp_global_def.print_timestamp("Finish")

if __name__ == "__main__":
    main()
