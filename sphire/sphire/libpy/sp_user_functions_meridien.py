
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
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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

#  This file contains fuctions that perform project-dependent tasks in various
#   alignment programs, for example preparation of the reference during 2D and 3D alignment
#  To write you own function, modify the existing one (for example, wei_func is a version
#   of ref_ali2d) and add the name to the factory.  Once it is done, the function can be called
#   from appropriate application, in this case "sxali2d_c.py ...  --function=wei_func
#

from . import sp_morphology
from . import sp_utilities
from . import sp_global_def
import past
from past.utils import old_div
import numpy
import EMAN2_cppwrap


def do_volume_mask(ref_data):
    """
		1. - volume
		2. - Tracker, see meridien
		3. - current iteration number
	"""

    # Retrieve the function specific input arguments from ref_data
    vol = ref_data[0]
    Tracker = ref_data[1]
    mainiteration = ref_data[2]

    if Tracker["constants"]["mask3D"] is None:
        vol = sp_morphology.cosinemask(vol, radius=Tracker["constants"]["radius"])
    else:
        EMAN2_cppwrap.Util.mul_img(
            vol, sp_utilities.get_im(Tracker["constants"]["mask3D"])
        )

    return vol


def compute_search_params(acc_trans, shifter, old_range):
    # step refer to the fine sampled step; while the range remains
    if old_range == 0.0 and shifter != 0.0:
        old_range = acc_trans
    step = min(1.5, 0.75 * acc_trans)  # remove 2
    new_range = min(1.3 * old_range, 5.0 * shifter)
    new_range = max(new_range, 3.0 * step)  # change 1.5 to 3.0
    if new_range > 8.0 * step:
        new_range = old_div(new_range, 2.0)  # change 4 to 8
    if new_range > 8.0 * step:
        step = old_div(new_range, 8.0)  # change 4 to 8
    #  change 4. to 8. on account of the fact that we return actual shift which is then doubled for coarse search.
    if new_range == 0.0:
        step = 0.5  # change 1.0 to 0.5
    return new_range, step


def ai_spa(Tracker, fff, anger, shifter, do_local, chout=False):
    """
	chout - if true, one can print, call the program with, chout = (Blockdata["myid"] == Blockdata["main_node"])
	fff (fsc), anger, shifter are coming from the previous iteration

	Possibilities we will consider:
	1.  resolution improved: keep going with current settings.
	2.  resolution stalled and no pwadjust: turn on pwadjust
	3.  resolution stalled and pwadjust: move to the next phase
	4.  resolution decreased: back off and move to the next phase
	5.  All phases tried and nxinit < nnxo: set nxinit == nnxo and run local searches.

	###  checkconvergence  merged in AI  02/16/2017
	# when the following conditions are all true
	#1. has_fine_enough_angular_sampling  True  #   Current sampling are fine enough
	#2. nr_iter_wo_resol_gain >= MAX_NR_ITER_WO_RESOL_GAIN #
	#3. nr_iter_wo_large_hidden_variable_changes >= MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES
	"""

    keepgoing = 1
    Tracker["keepfirst"] = -1

    l05 = -1
    l01 = -1
    for i in range(len(fff)):
        if fff[i] < 0.5:
            l05 = i - 1
            break
    for i in range(l05 + 1, len(fff)):
        if fff[i] < 0.143:
            l01 = i - 1
            break
    l01 = max(l01, -1)
    if fff:
        ai_string = (
            "  AI: Tracker[nxstep], TR[currentres], Tracker[fsc143], l05, l01, fff[old_div(Tracker[nxinit], 2)-1]:",
            Tracker["nxstep"],
            Tracker["currentres"],
            Tracker["fsc143"],
            l05,
            l01,
            fff[old_div(Tracker["nxinit"], 2) - 1],
        )

    if Tracker["mainiteration"] == 1 and not do_local:
        Tracker["state"] = "INITIAL"

        inc = Tracker["currentres"]
        if Tracker["large_at_Nyquist"]:
            inc += int(0.25 * past.utils.old_div(Tracker["constants"]["nnxo"], 2) + 0.5)
        else:
            inc += Tracker["nxstep"]
        Tracker["nxinit"] = int(
            min(2 * inc, Tracker["constants"]["nnxo"])
        )  #  Cannot exceed image size
        Tracker["local"] = False
        Tracker["changed_delta"] = False

    elif Tracker["mainiteration"] == 1 and do_local:
        if chout:
            sp_global_def.sxprint(ai_string)
        Tracker["state"] = "PRIMARY LOCAL"
        Tracker["currentres"] = l05
        Tracker["fsc143"] = l01
        Tracker["large_at_Nyquist"] = bool(
            fff[old_div(Tracker["nxinit"], 2)] > 0.1
            or fff[old_div(Tracker["nxinit"], 2) - 1] > 0.2
        )
        Tracker["nxinit"] = min(
            2 * Tracker["fsc143"], Tracker["constants"]["nnxo"]
        )  #  Cannot exceed image size
        Tracker["local"] = True
        Tracker["an"] = 6 * Tracker["delta"]
        Tracker["no_improvement"] = 0
        Tracker["no_params_changes"] = 0
        Tracker["anger"] = 1.0e23
        Tracker["shifter"] = 1.0e23
        Tracker["constants"]["best"] = Tracker["mainiteration"]
        if chout:
            sp_global_def.sxprint(
                "ITERATION  #%2d. Resolution achieved       : %3d/%3d pixels, %5.2fA/%5.2fA."
                % (
                    Tracker["mainiteration"],
                    Tracker["currentres"],
                    Tracker["fsc143"],
                    old_div(
                        Tracker["constants"]["pixel_size"]
                        * Tracker["constants"]["nnxo"],
                        float(Tracker["currentres"]),
                    ),
                    old_div(
                        Tracker["constants"]["pixel_size"]
                        * Tracker["constants"]["nnxo"],
                        float(Tracker["fsc143"]),
                    ),
                )
            )

    else:
        if chout:
            sp_global_def.sxprint(ai_string)

        if Tracker["mainiteration"] == 2 and not do_local:
            Tracker["state"] = "PRIMARY"

        if Tracker["mainiteration"] > 3 or not do_local:
            Tracker["nxstep"] = max(Tracker["nxstep"], l01 - l05 + 5)

        if Tracker["state"] == "FINAL" or Tracker["state"] == "RESTRICTED":
            Tracker["large_at_Nyquist"] = bool(
                fff[old_div(Tracker["nxinit"], 2)] > 0.1
                or fff[old_div(Tracker["nxinit"], 2) - 1] > 0.2
            )
        else:
            Tracker["large_at_Nyquist"] = bool(
                fff[old_div(Tracker["nxinit"], 2) - 1] > 0.2
            )

        if Tracker["mainiteration"] == 2 and not do_local:
            maxres = Tracker["constants"]["inires"]
            maxres_143 = l01
        else:
            maxres = max(
                l05, 5
            )  #  5 is minimum resolution of the map, could be set by the user
            maxres_143 = l01

        try:
            bestres_143 = Tracker["bestres_143"]
        except:
            Tracker["bestres_143"] = maxres_143

        if maxres >= Tracker["bestres"] and maxres_143 >= Tracker["bestres_143"]:
            Tracker["bestres"] = maxres
            Tracker["bestres_143"] = maxres_143
            Tracker["constants"]["best"] = max(Tracker["mainiteration"] - 1, 3)  #

        if maxres > Tracker["currentres"]:
            Tracker["no_improvement"] = 0
            Tracker["no_params_changes"] = 0
        else:
            Tracker["no_improvement"] += 1

        Tracker["currentres"] = maxres
        Tracker["fsc143"] = maxres_143

        params_changes = bool(
            anger >= 1.03 * Tracker["anger"] and shifter >= 1.03 * Tracker["shifter"]
        )

        #  figure changes in params
        if chout:
            sp_global_def.sxprint(
                "  Incoming  parameters  {0:10.3f}  {1:10.3f}  {2:10.3f}  {3:10.3f}   {4}".format(
                    Tracker["anger"], anger, Tracker["shifter"], shifter, params_changes
                )
            )
        if params_changes:
            Tracker["no_params_changes"] += 1
        else:
            Tracker["no_params_changes"] = 0

        if anger < Tracker["anger"]:
            Tracker["anger"] = anger
        if shifter < Tracker["shifter"]:
            Tracker["shifter"] = shifter

        inc = Tracker["currentres"]
        if Tracker["large_at_Nyquist"]:
            inc += int(0.25 * past.utils.old_div(Tracker["constants"]["nnxo"], 2) + 0.5)
            slim = int(Tracker["nxinit"] * 1.09)
            tmp = min(max(2 * inc, slim + slim % 2), Tracker["constants"]["nnxo"])
        else:
            inc += Tracker["nxstep"]
            tmp = min(
                2 * inc, Tracker["constants"]["nnxo"]
            )  #  Cannot exceed image size

        if chout:
            sp_global_def.sxprint(
                "  IN AI: nxstep, large at Nyq, outcoming current res, adjusted current, inc, estimated image size",
                Tracker["nxstep"],
                Tracker["large_at_Nyquist"],
                Tracker["currentres"],
                inc,
                tmp,
            )

        if do_local:
            tmp = max(tmp, Tracker["nxinit"])
        Tracker["nxinit"] = int(tmp)
        Tracker["changed_delta"] = False
        #  decide angular step and translations
        if (
            Tracker["no_improvement"] >= Tracker["constants"]["limit_improvement"]
            and Tracker["no_params_changes"] >= Tracker["constants"]["limit_changes"]
            and not Tracker["large_at_Nyquist"]
        ):
            if (
                Tracker["delta"]
                < Tracker["constants"]["a_criterion"] * Tracker["acc_rot"]
            ):  # <<<----it might cause converge issues when shake is 0.0
                if Tracker["state"] == "PRIMARY LOCAL":
                    step_range, step = compute_search_params(
                        Tracker["acc_trans"], Tracker["shifter"], Tracker["xr"]
                    )
                    if do_local:
                        step_range = min(step_range, Tracker["xr"])
                        step = min(step, Tracker["ts"])
                    if chout:
                        sp_global_def.sxprint(
                            "  Computed  pares  ",
                            Tracker["anger"],
                            anger,
                            Tracker["shifter"],
                            shifter,
                            Tracker["xr"],
                            step_range,
                            step,
                        )
                    Tracker["xr"] = step_range
                    Tracker["ts"] = step
                    Tracker["state"] = "RESTRICTED"
                else:
                    keepgoing = 0
                    if chout:
                        sp_global_def.sxprint(
                            "Convergence criterion A is reached (angular step delta smaller than 3/4 changes in angles))"
                        )
            else:
                step_range, step = compute_search_params(
                    Tracker["acc_trans"], Tracker["shifter"], Tracker["xr"]
                )
                if do_local:
                    step_range = min(step_range, Tracker["xr"])
                    step = min(step, Tracker["ts"])
                if chout:
                    sp_global_def.sxprint(
                        "  Computed  pares  ",
                        Tracker["anger"],
                        anger,
                        Tracker["shifter"],
                        shifter,
                        Tracker["xr"],
                        step_range,
                        step,
                    )
                Tracker["xr"] = step_range
                Tracker["ts"] = step
                Tracker["delta"] = old_div(Tracker["delta"], 2.0)
                Tracker["changed_delta"] = True
                if (
                    Tracker["delta"] <= old_div(3.75, 2.0) or do_local
                ):  #  MOVE DOWN TO RESTRICTED
                    Tracker["an"] = 6 * Tracker["delta"]
                    if Tracker["delta"] <= numpy.degrees(
                        numpy.arctan(old_div(0.25, Tracker["constants"]["radius"]))
                    ):
                        Tracker["state"] = "FINAL"
                    else:
                        Tracker["state"] = "RESTRICTED"
                else:
                    Tracker["an"] = -1
                    if Tracker["state"] == "PRIMARY":
                        Tracker["state"] = "EXHAUSTIVE"
                if chout:
                    sp_global_def.sxprint(
                        "  IN AI there was reset due to no changes, adjust stuff  ",
                        Tracker["no_improvement"],
                        Tracker["no_params_changes"],
                        Tracker["delta"],
                        Tracker["xr"],
                        Tracker["ts"],
                        Tracker["state"],
                    )
                # check convergence before reset
                if (
                    Tracker["state"] == "FINAL"
                    and Tracker["no_improvement"]
                    >= Tracker["constants"]["limit_improvement"]
                ):
                    keepgoing = 0
                    if chout:
                        sp_global_def.sxprint(
                            "Convergence criterion B is reached (angular step delta smaller than the limit imposed by the structure radius)"
                        )
                Tracker["no_improvement"] = 0
                Tracker["no_params_changes"] = 0
                Tracker["anger"] = 1.0e23
                Tracker["shifter"] = 1.0e23
    return keepgoing


def ai_filament(Tracker, fff, anger, shifter, do_local, chout=False):
    """
	chout - if true, one can print, call the program with, chout = (Blockdata["myid"] == Blockdata["main_node"])
	fff (fsc), anger, shifter are coming from the previous iteration

	Possibilities we will consider:
	1.  resolution improved: keep going with current settings.
	2.  resolution stalled and no pwadjust: turn on pwadjust
	3.  resolution stalled and pwadjust: move to the next phase
	4.  resolution decreased: back off and move to the next phase
	5.  All phases tried and nxinit < nnxo: set nxinit == nnxo and run local searches.

	###  checkconvergence  merged in AI  02/16/2017
	# when the following conditions are all true
	#1. has_fine_enough_angular_sampling  True  #   Current sampling are fine enough
	#2. nr_iter_wo_resol_gain >= MAX_NR_ITER_WO_RESOL_GAIN #
	#3. nr_iter_wo_large_hidden_variable_changes >= MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES
	"""

    keepgoing = 1
    Tracker["keepfirst"] = -1

    l05 = -1
    l01 = -1
    for i in range(len(fff)):
        if fff[i] < 0.5:
            l05 = i - 1
            break
    for i in range(l05 + 1, len(fff)):
        if fff[i] < 0.143:
            l01 = i - 1
            break
    l01 = max(l01, -1)
    if fff:
        ai_string = (
            "  AI: Tracker[nxstep], TR[currentres], Tracker[fsc143], l05, l01, fff[old_div(Tracker[nxinit], 2)-1]:",
            Tracker["nxstep"],
            Tracker["currentres"],
            Tracker["fsc143"],
            l05,
            l01,
            fff[old_div(Tracker["nxinit"], 2) - 1],
        )

    if Tracker["mainiteration"] == 1 and not do_local:
        Tracker["state"] = "INITIAL"

        inc = Tracker["currentres"]
        if Tracker["large_at_Nyquist"]:
            inc += int(0.25 * past.utils.old_div(Tracker["constants"]["nnxo"], 2) + 0.5)
        else:
            inc += Tracker["nxstep"]
        Tracker["nxinit"] = int(
            min(2 * inc, Tracker["constants"]["nnxo"])
        )  #  Cannot exceed image size
        Tracker["local"] = False
        Tracker["changed_delta"] = False

    elif Tracker["mainiteration"] == 1 and do_local:
        if chout:
            sp_global_def.sxprint(ai_string)
        Tracker["state"] = "PRIMARY LOCAL"
        Tracker["currentres"] = l05
        Tracker["fsc143"] = l01
        Tracker["large_at_Nyquist"] = bool(
            fff[old_div(Tracker["nxinit"], 2)] > 0.1
            or fff[old_div(Tracker["nxinit"], 2) - 1] > 0.2
        )
        Tracker["nxinit"] = min(
            2 * Tracker["fsc143"], Tracker["constants"]["nnxo"]
        )  #  Cannot exceed image size
        Tracker["local"] = True
        Tracker["an"] = 6 * Tracker["delta"]
        Tracker["no_improvement"] = 0
        Tracker["no_params_changes"] = 0
        Tracker["anger"] = 1.0e23
        Tracker["shifter"] = 1.0e23
        Tracker["constants"]["best"] = Tracker["mainiteration"]
        if chout:
            sp_global_def.sxprint(
                "ITERATION  #%2d. Resolution achieved       : %3d/%3d pixels, %5.2fA/%5.2fA."
                % (
                    Tracker["mainiteration"],
                    Tracker["currentres"],
                    Tracker["fsc143"],
                    old_div(
                        Tracker["constants"]["pixel_size"]
                        * Tracker["constants"]["nnxo"],
                        float(Tracker["currentres"]),
                    ),
                    old_div(
                        Tracker["constants"]["pixel_size"]
                        * Tracker["constants"]["nnxo"],
                        float(Tracker["fsc143"]),
                    ),
                )
            )

    else:
        if chout:
            sp_global_def.sxprint(ai_string)

        if Tracker["mainiteration"] == 2 and not do_local:
            Tracker["state"] = "PRIMARY"
        elif (
            Tracker["mainiteration"] == 6
            and not do_local
            and Tracker["state"] == "PRIMARY"
        ):
            Tracker["state"] = "EXHAUSTIVE"
        elif (
            Tracker["mainiteration"] == 12
            and not do_local
            and Tracker["state"] == "EXHAUSTIVE"
            and Tracker["delta"] <= 3.75
        ):
            Tracker["state"] = "RESTRICTED"
            Tracker["an"] = 6 * Tracker["delta"]
            Tracker["theta_min"] = 40
            Tracker["theta_max"] = 140
            Tracker["constants"]["shake"] = 0.5
            Tracker["delta"] = old_div(Tracker["delta"], 2.0)

        if Tracker["mainiteration"] > 3 or not do_local:
            Tracker["nxstep"] = max(Tracker["nxstep"], l01 - l05 + 5)

        if Tracker["state"] == "FINAL" or Tracker["state"] == "RESTRICTED":
            Tracker["large_at_Nyquist"] = bool(
                fff[old_div(Tracker["nxinit"], 2)] > 0.1
                or fff[old_div(Tracker["nxinit"], 2) - 1] > 0.2
            )
        else:
            Tracker["large_at_Nyquist"] = bool(
                fff[old_div(Tracker["nxinit"], 2) - 1] > 0.2
            )

        if Tracker["mainiteration"] == 2 and not do_local:
            maxres = Tracker["constants"]["inires"]
            maxres_143 = l01
        else:
            maxres = max(
                l05, 5
            )  #  5 is minimum resolution of the map, could be set by the user
            maxres_143 = l01

        try:
            bestres_143 = Tracker["bestres_143"]
        except:
            Tracker["bestres_143"] = maxres_143

        if maxres >= Tracker["bestres"] and maxres_143 >= Tracker["bestres_143"]:
            Tracker["bestres"] = maxres
            Tracker["bestres_143"] = maxres_143
            Tracker["constants"]["best"] = max(Tracker["mainiteration"] - 1, 3)  #

        if maxres > Tracker["currentres"]:
            Tracker["no_improvement"] = 0
            Tracker["no_params_changes"] = 0
        else:
            Tracker["no_improvement"] += 1

        Tracker["currentres"] = maxres
        Tracker["fsc143"] = maxres_143

        params_changes = bool(
            anger >= 1.03 * Tracker["anger"] and shifter >= 1.03 * Tracker["shifter"]
        )

        #  figure changes in params
        if chout:
            sp_global_def.sxprint(
                "  Incoming  parameters  {0:10.3f}  {1:10.3f}  {2:10.3f}  {3:10.3f}   {4}".format(
                    Tracker["anger"], anger, Tracker["shifter"], shifter, params_changes
                )
            )
        if params_changes:
            Tracker["no_params_changes"] += 1
        else:
            Tracker["no_params_changes"] = 0

        if anger < Tracker["anger"]:
            Tracker["anger"] = anger
        if shifter < Tracker["shifter"]:
            Tracker["shifter"] = shifter

        inc = Tracker["currentres"]
        if Tracker["large_at_Nyquist"]:
            inc += int(0.25 * past.utils.old_div(Tracker["constants"]["nnxo"], 2) + 0.5)
            slim = int(Tracker["nxinit"] * 1.09)
            tmp = min(max(2 * inc, slim + slim % 2), Tracker["constants"]["nnxo"])
        else:
            inc += Tracker["nxstep"]
            tmp = min(
                2 * inc, Tracker["constants"]["nnxo"]
            )  #  Cannot exceed image size

        if chout:
            sp_global_def.sxprint(
                "  IN AI: nxstep, large at Nyq, outcoming current res, adjusted current, inc, estimated image size",
                Tracker["nxstep"],
                Tracker["large_at_Nyquist"],
                Tracker["currentres"],
                inc,
                tmp,
            )

        tmp = max(tmp, Tracker["nxinit"])
        Tracker["nxinit"] = int(tmp)
        Tracker["changed_delta"] = False
        #  decide angular step and translations
        if (
            Tracker["no_improvement"] >= Tracker["constants"]["limit_improvement"]
            and Tracker["no_params_changes"] >= Tracker["constants"]["limit_changes"]
            and not Tracker["large_at_Nyquist"]
        ):
            if (
                Tracker["delta"]
                < Tracker["constants"]["a_criterion"] * Tracker["acc_rot"]
            ):  # <<<----it might cause converge issues when shake is 0.0
                if Tracker["state"] == "PRIMARY LOCAL":
                    step_range, step = compute_search_params(
                        Tracker["acc_trans"], Tracker["shifter"], Tracker["xr"]
                    )
                    if do_local:
                        step_range = min(step_range, Tracker["xr"])
                        step = min(step, Tracker["ts"])
                    if chout:
                        sp_global_def.sxprint(
                            "  Computed  pares  ",
                            Tracker["anger"],
                            anger,
                            Tracker["shifter"],
                            shifter,
                            Tracker["xr"],
                            step_range,
                            step,
                        )
                    Tracker["xr"] = step_range
                    Tracker["ts"] = step
                    Tracker["state"] = "RESTRICTED"
                else:
                    keepgoing = 0
                    if chout:
                        sp_global_def.sxprint(
                            "Convergence criterion A is reached (angular step delta smaller than 3/4 changes in angles))"
                        )
            else:
                step_range, step = compute_search_params(
                    Tracker["acc_trans"], Tracker["shifter"], Tracker["xr"]
                )
                step_range = min(step_range, Tracker["xr"])
                step = min(step, Tracker["ts"])
                if chout:
                    sp_global_def.sxprint(
                        "  Computed  pares  ",
                        Tracker["anger"],
                        anger,
                        Tracker["shifter"],
                        shifter,
                        Tracker["xr"],
                        step_range,
                        step,
                    )
                Tracker["xr"] = step_range
                Tracker["ts"] = step
                Tracker["delta"] = old_div(Tracker["delta"], 2.0)
                Tracker["changed_delta"] = True
                if Tracker["state"] == "PRIMARY":
                    Tracker["delta"] *= 2.0
                    Tracker["state"] = "EXHAUSTIVE"
                elif (
                    Tracker["delta"] <= old_div(3.75, 2.0) or do_local
                ):  #  MOVE DOWN TO RESTRICTED
                    Tracker["an"] = 6 * Tracker["delta"]
                    Tracker["theta_min"] = 40
                    Tracker["theta_max"] = 140
                    Tracker["constants"]["shake"] = 0.5
                    if Tracker["delta"] <= numpy.degrees(
                        numpy.arctan(old_div(0.25, Tracker["constants"]["radius"]))
                    ):
                        Tracker["state"] = "FINAL"
                    else:
                        Tracker["state"] = "RESTRICTED"
                else:
                    Tracker["an"] = -1
                    if Tracker["state"] == "PRIMARY":
                        Tracker["state"] = "EXHAUSTIVE"
                if chout:
                    sp_global_def.sxprint(
                        "  IN AI there was reset due to no changes, adjust stuff  ",
                        Tracker["no_improvement"],
                        Tracker["no_params_changes"],
                        Tracker["delta"],
                        Tracker["xr"],
                        Tracker["ts"],
                        Tracker["state"],
                    )
                # check convergence before reset
                if (
                    Tracker["state"] == "FINAL"
                    and Tracker["no_improvement"]
                    >= Tracker["constants"]["limit_improvement"]
                ):
                    keepgoing = 0
                    if chout:
                        sp_global_def.sxprint(
                            "Convergence criterion B is reached (angular step delta smaller than the limit imposed by the structure radius)"
                        )
                Tracker["no_improvement"] = 0
                Tracker["no_params_changes"] = 0
                Tracker["anger"] = 1.0e23
                Tracker["shifter"] = 1.0e23

    if Tracker["state"] == "RESTRICTED":
        Tracker["constants"]["do_rotate"] = True
        Tracker["ccfpercentage"] = min(Tracker["ccfpercentage"] + 0.2, 0.999)
        Tracker["prior"]["force_outlier"] = False
        Tracker["prior"]["apply_prior"] = True
    elif Tracker["state"] == "EXHAUSTIVE":
        Tracker["ccfpercentage"] = 0.3
        Tracker["constants"]["do_rotate"] = False
        Tracker["prior"]["force_outlier"] = False
        Tracker["prior"]["apply_prior"] = True
    else:
        Tracker["constants"]["do_rotate"] = False
        Tracker["prior"]["force_outlier"] = False
        Tracker["prior"]["apply_prior"] = True

    return keepgoing
