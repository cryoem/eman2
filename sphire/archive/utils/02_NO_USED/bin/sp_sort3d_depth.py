



















































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































from __future__ import print_function
def get_shrink_data_sorting_smearing(partids, partstack, return_real=False, preshift=True, apply_mask=True, npad=1):
    # The function will read from stack a subset of images specified in partids
    #   and assign to them parameters from partstack with optional CTF application and shifting of the data.
    # So, the lengths of partids and partstack are the same.
    #  The read data is properly distributed among MPI threads.
    # 10142015 --- preshift is set to True when doing 3-D sorting.
    # chunk_id are set when data is read in
    global Tracker, Blockdata
    pass#IMPORTIMPORTIMPORT from sp_fundamentals import resample, fshift
    pass#IMPORTIMPORTIMPORT from sp_filter import filt_ctf
    pass#IMPORTIMPORTIMPORT from sp_applications import MPI_start_end
    pass#IMPORTIMPORTIMPORT from EMAN2 import Region
    mask2D = sp_utilities.model_circle(Tracker["constants"]["radius"], Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"])
    shrinkage = Tracker["nxinit"] / float(Tracker["constants"]["nnxo"])
    radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)

    if (Blockdata["myid"] == Blockdata["main_node"]):
        lpartids = sp_utilities.read_text_file(partids, -1)
        if len(lpartids) == 1:
            lpartids = lpartids[0]
            groupids = len(lpartids) * [-1]
        else:
            groupids = lpartids[0]
            lpartids = lpartids[1]
    else:
        lpartids = 0
        groupids = 0
    lpartids = sp_utilities.wrap_mpi_bcast(lpartids, Blockdata["main_node"])
    groupids = sp_utilities.wrap_mpi_bcast(groupids, Blockdata["main_node"])
    Tracker["total_stack"] = len(lpartids)
    if (Blockdata["myid"] == Blockdata["main_node"]):
        partstack = sp_utilities.read_text_row(partstack)
    else:
        partstack = 0
    partstack = sp_utilities.wrap_mpi_bcast(partstack, Blockdata["main_node"])
    if (Tracker["total_stack"] < Blockdata["nproc"]):
        sp_global_def.ERROR("Wrong MPI settings!", "get_shrink_data_sorting", 1, Blockdata["myid"])
    else:
        image_start, image_end = sp_applications.MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
    lpartids = lpartids[image_start:image_end]
    groupids = groupids[image_start:image_end]
    nima = image_end - image_start
    data = [None] * nima
    norm_per_particle = []
    for im in range(nima):
        data[im] = sp_utilities.get_im(Tracker["constants"]["orgstack"], lpartids[im])
        try:
            phi, theta, psi, sx, sy, chunk_id, particle_group_id, norm = partstack[lpartids[im]][0], \
                                                                         partstack[lpartids[im]][1], \
                                                                         partstack[lpartids[im]][2], \
                                                                         partstack[lpartids[im]][3], \
                                                                         partstack[lpartids[im]][4], \
                                                                         partstack[lpartids[im]][5], \
                                                                         partstack[lpartids[im]][6], \
                                                                         partstack[lpartids[im]][7]
        except:
            phi, theta, psi, sx, sy, chunk_id, particle_group_id, norm = partstack[lpartids[im]][0], \
                                                                         partstack[lpartids[im]][1], \
                                                                         partstack[lpartids[im]][2], \
                                                                         partstack[lpartids[im]][3], \
                                                                         partstack[lpartids[im]][4], \
                                                                         partstack[lpartids[im]][5], -1, 1
        if preshift:  # always true
            data[im] = sp_fundamentals.fshift(data[im], sx, sy)
            sx = 0.0
            sy = 0.0
        st = EMAN2_cppwrap.Util.infomask(data[im], mask2D, False)
        data[im] -= st[0]
        data[im] /= st[1]
        if apply_mask: data[im] = sp_morphology.cosinemask(data[im], radius=Tracker["constants"]["radius"])
        # FT
        data[im] = sp_fundamentals.fft(data[im])
        nny = data[im].get_ysize()
        if Tracker["constants"]["CTF"]:
            ctf_params = data[im].get_attr("ctf")
            data[im] = sp_fundamentals.fdecimate(data[im], Tracker["nxinit"] * npad, Tracker["nxinit"] * npad, 1, False, False)
            ctf_params.apix = ctf_params.apix / shrinkage
            data[im].set_attr('ctf', ctf_params)
            data[im].set_attr('ctf_applied', 0)
            if return_real:  data[im] = sp_fundamentals.fft(data[im])
        else:
            ctf_params = data[im].get_attr_default("ctf", False)
            if ctf_params:
                ctf_params.apix = ctf_params.apix / shrinkage
                data[im].set_attr('ctf', ctf_params)
                data[im].set_attr('ctf_applied', 0)
            data[im] = sp_fundamentals.fdecimate(data[im], nxinit * npad, nxinit * npad, 1, True, False)
            apix = Tracker["constants"]["pixel_size"]
            data[im].set_attr('apix', apix / shrinkage)
        if not return_real:    data[im].set_attr("padffted", 1)
        data[im].set_attr("npad", npad)
        sp_utilities.set_params_proj(data[im], [phi, theta, psi, 0.0, 0.0])
        data[im].set_attr("chunk_id", chunk_id)
        data[im].set_attr("group", groupids[im])
        data[im].set_attr("particle_group", particle_group_id)
        if Tracker["applybckgnoise"]:
            data[im].set_attr("bckgnoise", Blockdata["bckgnoise"][particle_group_id])
            data[im].set_attr("qt", float(Tracker["constants"]["nnxo"] * Tracker["constants"]["nnxo"]))
        else:
            data[im].set_attr("bckgnoise", Blockdata["bckgnoise"])  # constant list
        norm_per_particle.append(norm)
    return data, norm_per_particle


###3
def get_data_prep_compare_rec3d(partids, partstack, return_real=False, preshift=True, npad=1):
    # The function will read from stack a subset of images specified in partids
    #   and assign to them parameters from partstack with optional CTF application and shifting of the data.
    # So, the lengths of partids and partstack are the same.

    global Tracker, Blockdata
    pass#IMPORTIMPORTIMPORT from sp_fundamentals import resample, fshift, fft
    pass#IMPORTIMPORTIMPORT from sp_filter import filt_ctf
    pass#IMPORTIMPORTIMPORT from sp_applications import MPI_start_end
    pass#IMPORTIMPORTIMPORT from EMAN2 import Region
    pass#IMPORTIMPORTIMPORT from sp_utilities import model_circle, wrap_mpi_bcast, get_im, model_blank, set_params_proj
    # functions:
    # read in data
    # apply mask, and prepare focus projection if focus3D is specified
    # return  1. cdata: data for image comparison, always in Fourier format
    #         2. rdata: data for reconstruction, 4nn return real image

    mask2D = sp_utilities.model_circle(Tracker["constants"]["radius"], Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"])
    shrinkage = Tracker["nxinit"] / float(Tracker["constants"]["nnxo"])
    radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
    if Tracker["applybckgnoise"]:
        oneover = []
        nnx = len(Blockdata["bckgnoise"][0])
        for i in range(len(Blockdata["bckgnoise"])):
            temp = [0.0] * nnx
            for k in range(nnx):
                if (Blockdata["bckgnoise"][i][k] > 0.0):  temp[k] = 1.0 / numpy.sqrt(Blockdata["bckgnoise"][i][k])
            oneover.append(temp)
        del temp
    if (Blockdata["myid"] == Blockdata["main_node"]):
        lpartids = sp_utilities.read_text_file(partids, -1)
        if len(lpartids) == 1:
            lpartids = lpartids[0]
            groupids = len(lpartids) * [-1]
        else:
            groupids = lpartids[0]
            lpartids = lpartids[1]
    else:
        lpartids = 0
        groupids = 0
    lpartids = sp_utilities.wrap_mpi_bcast(lpartids, Blockdata["main_node"])
    groupids = sp_utilities.wrap_mpi_bcast(groupids, Blockdata["main_node"])
    Tracker["total_stack"] = len(lpartids)
    if (Blockdata["myid"] == Blockdata["main_node"]):
        partstack = sp_utilities.read_text_row(partstack)
    else:
        partstack = 0
    partstack = sp_utilities.wrap_mpi_bcast(partstack, Blockdata["main_node"])
    if (Tracker["total_stack"] < Blockdata["nproc"]):
        sp_global_def.ERROR("Number of processors is larger than the total number of images", \
              "get_data_and_prep", 1, Blockdata["myid"])
    else:
        image_start, image_end = sp_applications.MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])
    lpartids = lpartids[image_start:image_end]
    groupids = groupids[image_start:image_end]
    if Tracker["constants"]["focus3D"]:  # focus mask is applied
        if Blockdata["myid"] == Blockdata["main_node"]:
            focus3d = sp_utilities.get_im(Tracker["constants"]["focus3D"])
            focus3d_nx = focus3d.get_xsize()
            if focus3d_nx != Tracker["constants"]["nnxo"]:  # So the decimated focus volume can be directly used
                focus3d = sp_fundamentals.resample(focus3d, float(Tracker["constants"]["nnxo"]) / float(focus3d_nx))
        else:
            focus3d = sp_utilities.model_blank(Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"],
                                  Tracker["constants"]["nnxo"])
        sp_utilities.bcast_EMData_to_all(focus3d, Blockdata["myid"], Blockdata["main_node"])
        focus3d = sp_projection.prep_vol(focus3d, 1, 1)
    #  Preprocess the data
    #  mask2D    =    model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
    nima = image_end - image_start
    cdata = [None] * nima
    rdata = [None] * nima
    for im in range(nima):
        image = sp_utilities.get_im(Tracker["constants"]["orgstack"], lpartids[im])
        try:
            phi, theta, psi, sx, sy, chunk_id, particle_group_id = partstack[lpartids[im]][0], partstack[lpartids[im]][
                1], partstack[lpartids[im]][2], \
                                                                   partstack[lpartids[im]][3], partstack[lpartids[im]][
                                                                       4], partstack[lpartids[im]][5], \
                                                                   partstack[lpartids[im]][6]
        except:
            phi, theta, psi, sx, sy, chunk_id, particle_group_id = partstack[lpartids[im]][0], partstack[lpartids[im]][
                1], partstack[lpartids[im]][2], \
                                                                   partstack[lpartids[im]][3], partstack[lpartids[im]][
                                                                       4], partstack[lpartids[im]][5], -1
        if preshift:  # always true
            image = sp_fundamentals.fshift(image, sx, sy)
            sx = 0.0
            sy = 0.0
        st = EMAN2_cppwrap.Util.infomask(image, mask2D, False)
        image -= st[0]
        image /= st[1]
        cimage = image.copy()
        if Tracker["applybckgnoise"]:
            if Tracker["applymask"]:
                if Tracker["constants"]["hardmask"]:
                    cimage = sp_morphology.cosinemask(cimage, radius=Tracker["constants"]["radius"])
                else:
                    bckg = sp_utilities.model_gauss_noise(1.0, Tracker["constants"]["nnxo"] + 2, Tracker["constants"]["nnxo"])
                    bckg.set_attr("is_complex", 1)
                    bckg.set_attr("is_fftpad", 1)
                    bckg = sp_fundamentals.fft(sp_filter.filt_table(bckg, oneover[particle_group_id]))
                    #  Normalize bckg noise in real space, only region actually used.
                    st = EMAN2_cppwrap.Util.infomask(bckg, mask2D, False)
                    bckg -= st[0]
                    bckg /= st[1]
                    cimage = sp_morphology.cosinemask(cimage, radius=Tracker["constants"]["radius"], bckg=bckg)
        else:
            if Tracker["applymask"]: cimage = sp_morphology.cosinemask(cimage, radius=Tracker["constants"]["radius"])
        # FT
        image = sp_fundamentals.fft(image)
        cimage = sp_fundamentals.fft(cimage)
        if Tracker["constants"]["CTF"]:
            ctf_params = image.get_attr("ctf")
            image = sp_fundamentals.fdecimate(image, Tracker["nxinit"] * npad, Tracker["nxinit"] * npad, 1, False, False)
            cimage = sp_fundamentals.fdecimate(cimage, Tracker["nxinit"] * npad, Tracker["nxinit"] * npad, 1, False, False)
            ctf_params.apix = ctf_params.apix / shrinkage
            image.set_attr('ctf', ctf_params)
            cimage.set_attr('ctf', ctf_params)
            image.set_attr('ctf_applied', 0)
            cimage.set_attr('ctf_applied', 0)
            if return_real: image = sp_fundamentals.fft(image)
        else:
            ctf_params = image.get_attr_default("ctf", False)
            if ctf_params:
                ctf_params.apix = ctf_params.apix / shrinkage
                image.set_attr('ctf', ctf_params)
                image.set_attr('ctf_applied', 0)
                cimage.set_attr('ctf', ctf_params)
                cimage.set_attr('ctf_applied', 0)
            image = sp_fundamentals.fdecimate(image, nxinit * npad, nxinit * npad, 1, True, False)
            cimage = sp_fundamentals.fdecimate(cimage, nxinit * npad, nxinit * npad, 1, True, False)
            apix = Tracker["constants"]["pixel_size"]
            image.set_attr('apix', apix / shrinkage)
            cimage.set_attr('apix', apix / shrinkage)
        cimage.set_attr("padffted", 1)
        cimage.set_attr("npad", npad)
        if not return_real:
            image.set_attr("padffted", 1)
            image.set_attr("npad", npad)
        sp_utilities.set_params_proj(image, [phi, theta, psi, 0.0, 0.0])
        image.set_attr("chunk_id", chunk_id)
        image.set_attr("group", groupids[im])
        image.set_attr("particle_group", particle_group_id)
        sp_utilities.set_params_proj(cimage, [phi, theta, psi, 0.0, 0.0])
        cimage.set_attr("chunk_id", chunk_id)
        cimage.set_attr("group", groupids[im])
        cimage.set_attr("particle_group", particle_group_id)
        rdata[im] = image
        cdata[im] = cimage
        if Tracker["applybckgnoise"]:
            rdata[im].set_attr("bckgnoise", Blockdata["bckgnoise"][particle_group_id])
            if Tracker["constants"]["comparison_method"] == "cross": EMAN2_cppwrap.Util.mulclreal(cdata[im], Blockdata["unrolldata"][
                particle_group_id])
        if Tracker["constants"]["focus3D"]:
            cdata[im] = sp_fundamentals.fft(sp_morphology.binarize(sp_projection.prgl(focus3d, [phi, theta, psi, 0.0, 0.0], 1, True), 1) * sp_fundamentals.fft(cdata[im]))
            if Tracker["constants"]["CTF"]: cdata[im].set_attr("ctf", rdata[im].get_attr("ctf"))
        cdata[im].set_attr("is_complex", 0)
    return cdata, rdata


#####4
def get_shrink_data_final(nxinit, procid, original_data=None, oldparams=None, \
                          return_real=False, preshift=False, apply_mask=True, nonorm=False, npad=1):
    global Tracker, Blockdata
    """
    This function will read from stack a subset of images specified in partids
       and assign to them parameters from partstack with optional CTF application and shifting of the data.
    So, the lengths of partids and partstack are the same.
      The read data is properly distributed among MPI threads.

    Flow of data:
    1. Read images, if there is enough memory, keep them as original_data.
    2. Read current params
    3.  Apply shift
    4.  Normalize outside of the radius
    5.  Do noise substitution and cosine mask.  (Optional?)
    6.  Shrink data.
    7.  Apply CTF.

    """
    # from fundamentals import resample
    pass#IMPORTIMPORTIMPORT from sp_utilities import get_im, model_gauss_noise, set_params_proj, get_params_proj
    pass#IMPORTIMPORTIMPORT from sp_fundamentals import fdecimate, fshift, fft
    pass#IMPORTIMPORTIMPORT from sp_filter import filt_ctf, filt_table
    pass#IMPORTIMPORTIMPORT from sp_applications import MPI_start_end
    pass#IMPORTIMPORTIMPORT from math import sqrt

    mask2D = sp_utilities.model_circle(Tracker["constants"]["radius"], Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"])
    nima = len(original_data)
    shrinkage = nxinit / float(Tracker["constants"]["nnxo"])
    radius = int(Tracker["constants"]["radius"] * shrinkage + 0.5)
    txm = float(nxinit - (nxinit // 2 + 1) - radius)
    txl = float(radius - nxinit // 2 + 1)

    if Blockdata["bckgnoise"]:
        oneover = []
        nnx = Blockdata["bckgnoise"][0].get_xsize()
        for i in range(len(Blockdata["bckgnoise"])):
            temp = [0.0] * nnx
            for k in range(nnx):
                if (Blockdata["bckgnoise"][i].get_value_at(k) > 0.0):  temp[k] = 1.0 / numpy.sqrt(
                    Blockdata["bckgnoise"][i].get_value_at(k))
            oneover.append(temp)
        del temp
    Blockdata["accumulatepw"][procid] = [None] * nima
    data = [None] * nima
    for im in range(nima):
        phi, theta, psi, sx, sy, wnorm = oldparams[im][0], oldparams[im][1], oldparams[im][2], oldparams[im][3], \
                                         oldparams[im][4], oldparams[im][7]
        if preshift:
            sx = int(round(sx))
            sy = int(round(sy))
            data[im] = sp_fundamentals.cyclic_shift(original_data[im], sx, sy)
            #  Put rounded shifts on the list, note image has the original floats - check whether it may cause problems
            oldparams[im][3] = sx
            oldparams[im][4] = sy
            sx = 0.0
            sy = 0.0
        else:
            data[im] = original_data[im].copy()
        st = EMAN2_cppwrap.Util.infomask(data[im], mask2D, False)
        data[im] -= st[0]
        data[im] /= st[1]
        if data[im].get_attr_default("bckgnoise", None):  data[im].delete_attr("bckgnoise")
        #  Do bckgnoise if exists
        if Blockdata["bckgnoise"]:
            if apply_mask:
                if Tracker["constants"]["hardmask"]:
                    data[im] = sp_morphology.cosinemask(data[im], radius=Tracker["constants"]["radius"])
                else:
                    bckg = sp_utilities.model_gauss_noise(1.0, Tracker["constants"]["nnxo"] + 2, Tracker["constants"]["nnxo"])
                    bckg.set_attr("is_complex", 1)
                    bckg.set_attr("is_fftpad", 1)
                    bckg = sp_fundamentals.fft(sp_filter.filt_table(bckg, oneover[data[im].get_attr("particle_group")]))
                    #  Normalize bckg noise in real space, only region actually used.
                    st = EMAN2_cppwrap.Util.infomask(bckg, mask2D, False)
                    bckg -= st[0]
                    bckg /= st[1]
                    data[im] = sp_morphology.cosinemask(data[im], radius=Tracker["constants"]["radius"], bckg=bckg)
        else:
            #  if no bckgnoise, do simple masking instead
            if apply_mask:  data[im] = sp_morphology.cosinemask(data[im], radius=Tracker["constants"]["radius"])
        #  Apply varadj
        if not nonorm: EMAN2_cppwrap.Util.mul_scalar(data[im], Tracker["avgvaradj"][procid] / wnorm)
        #  FT
        data[im] = sp_fundamentals.fft(data[im])
        sig = EMAN2_cppwrap.Util.rotavg_fourier(data[im])
        Blockdata["accumulatepw"][procid][im] = sig[len(sig) // 2:] + [0.0]
        if Tracker["constants"]["CTF"]:
            data[im] = sp_fundamentals.fdecimate(data[im], nxinit * npad, nxinit * npad, 1, False, False)
            ctf_params = original_data[im].get_attr("ctf")
            ctf_params.apix = ctf_params.apix / shrinkage
            data[im].set_attr('ctf', ctf_params)
            data[im].set_attr('ctf_applied', 0)
            if return_real: data[im] = sp_fundamentals.fft(data[im])
        else:
            ctf_params = original_data[im].get_attr_default("ctf", False)
            if ctf_params:
                ctf_params.apix = ctf_params.apix / shrinkage
                data[im].set_attr('ctf', ctf_params)
                data[im].set_attr('ctf_applied', 0)
            data[im] = sp_fundamentals.fdecimate(data[im], nxinit * npad, nxinit * npad, 1, True, False)
            apix = Tracker["constants"]["pixel_size"]
            data[im].set_attr('apix', apix / shrinkage)

        #  We have to make sure the shifts are within correct range, shrinkage or not
        sp_utilities.set_params_proj(data[im],
                        [phi, theta, psi, max(min(sx * shrinkage, txm), txl), max(min(sy * shrinkage, txm), txl)])
        if not return_real: data[im].set_attr("padffted", 1)
        data[im].set_attr("npad", npad)
        if Blockdata["bckgnoise"]:
            temp = Blockdata["bckgnoise"][data[im].get_attr("particle_group")]
            ###  Do not adjust the values, we try to keep everything in the same Fourier values.
            data[im].set_attr("bckgnoise", [temp[i] for i in range(temp.get_xsize())])
    return data


###5











































































































































































































































































































































































































































































































































































































































































































































def find_smallest_group(clusters):
    min_size = [len(clusters[0]), [0]]
    for ic in range(1, len(clusters)):
        if len(cluster[ic]) < min_size[0]:
            min_size = [len(clusters[ic]), [ic]]
        elif len(cluster[ic]) == min_size[0]:
            min_size[1].append(ic)
    if len(min_size[1]) >= 1: random.shuffle(min_size[1])
    return min_size[1][0]










































































































































































def even_assignment_alist_to_mclusters(glist, number_of_groups):
    # evenly assign glist to clusters
    pass#IMPORTIMPORTIMPORT import copy
    if number_of_groups > 0:
        clusters = [[] for i in range(number_of_groups)]
        ulist = copy.copy.deepcopy(glist)
        nc = 0
        while len(ulist) > 0:
            im = nc % number_of_groups
            random.shuffle(ulist)
            clusters[im].append(ulist[0])
            del ulist[0]
            nc += 1
        return clusters
    else:
        return []

























































































































































































































































































































































































































































































































































































































































































































def convertasi(asig, number_of_groups):
    pass#IMPORTIMPORTIMPORT from numpy import array
    p = []
    for k in range(number_of_groups):
        l = []
        for i in range(len(asig)):
            if (asig[i] == k): l.append(i)
        l = numpy.array(l, "int32")
        l.sort()
        p.append(l)
    return p


def partition_data_into_orientation_groups_nompi(refa_vecs, data_vecs):
    orien_assignment = [None for im in range(len(data_vecs))]
    for im in range(len(data_vecs)):
        max_dist = -999.0
        for jm in range(len(refa_vecs)):
            this_dis = get_dist1(data_vecs[im], refa_vecs[jm])
            if this_dis > max_dist:
                max_dist = this_dis
                orien_assignment[im] = jm
    return orien_assignment


### dmatrix and refangles partition
def get_dist1(vec1, vec2):
    sum_dot = 0.0
    for icomp in range(len(vec1)): sum_dot += vec1[icomp] * vec2[icomp]
    return sum_dot


def find_neighborhood(refa_vecs, minor_groups):
    matched_oriens = [[None, None] for i in range(len(minor_groups))]
    for iproj in range(len(minor_groups)):
        max_value = -999.0
        for jproj in range(len(refa_vecs)):
            if jproj not in minor_groups:
                this_dis = get_dist1(refa_vecs[minor_groups[iproj]], refa_vecs[jproj])
                if this_dis > max_value:
                    max_value = this_dis
                    matched_oriens[iproj] = [minor_groups[iproj], jproj]
    return matched_oriens


def reassign_ptls_in_orien_groups(assigned_ptls_in_groups, matched_pairs):
    tmplist = []
    for iorien in range(len(matched_pairs)):
        if matched_pairs[iorien][1] != None and matched_pairs[iorien][0] != None:
            assigned_ptls_in_groups[matched_pairs[iorien][1]] += assigned_ptls_in_groups[matched_pairs[iorien][0]]
            tmplist.append(matched_pairs[iorien][0])
    reassignment = []
    for iorien in range(len(assigned_ptls_in_groups)):
        if iorien not in tmplist: reassignment.append(sorted(assigned_ptls_in_groups[iorien]))
    return reassignment


def findall_dict(value, L, start=0):
    """
     return a list of all indices of a value on the list L beginning from position start
    """
    positions = []
    lL = len(L) - 1
    i = start - 1
    while (i < lL):
        i += 1
        try:
            if value == L[i]: positions.append(i)
        except:
            pass
    return positions


def get_orien_assignment_mpi(angle_step, partids, params, log_main):
    global Tracker, Blockdata
    pass#IMPORTIMPORTIMPORT from sp_applications import MPI_start_end
    pass#IMPORTIMPORTIMPORT from sp_utilities import wrap_mpi_recv, wrap_mpi_bcast, wrap_mpi_send, read_text_row, read_text_file, getvec
    sym_class = Blockdata["symclass"]
    image_start, image_end = sp_applications.MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], Blockdata["myid"])

    if Blockdata["myid"] == Blockdata["main_node"]:
        orien_group_assignment = [None for im in range(Tracker["total_stack"])]
    else:
        orien_group_assignment = 0
    refa = sym_class.even_angles(angle_step, theta1=Tracker["tilt1"], theta2=Tracker["tilt2"])

    refa_vecs = []
    for i in range(len(refa)): refa_vecs.append(sp_utilities.getvec(refa[i][0], refa[i][1]))

    if Blockdata["main_node"] == Blockdata["myid"]:
        params = sp_utilities.read_text_row(params)
        partids = sp_utilities.read_text_file(partids, -1)
        if len(partids) == 1:
            partids = partids[0]
        else:
            partids = partids[1]
        data_angles = [[None, None] for im in range(len(partids))]
        for im in range(len(partids)):
            data_angles[im] = sp_utilities.getvec(params[partids[im]][0], params[partids[im]][1])
        del params
        del partids
    else:
        data_angles = 0
    data_angles = sp_utilities.wrap_mpi_bcast(data_angles, Blockdata["main_node"], mpi.MPI_COMM_WORLD)

    data_angles = data_angles[image_start: image_end]
    local_orien_group_assignment = partition_data_into_orientation_groups_nompi(refa_vecs, data_angles)
    if Blockdata["myid"] == Blockdata["main_node"]:
        orien_group_assignment[image_start:image_end] = local_orien_group_assignment[:]
    else:
        orien_group_assignment = 0
    if Blockdata["main_node"] != Blockdata["myid"]:
        sp_utilities.wrap_mpi_send(local_orien_group_assignment, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
    else:
        for iproc in range(Blockdata["nproc"]):
            iproc_image_start, iproc_image_end = sp_applications.MPI_start_end(Tracker["total_stack"], Blockdata["nproc"], iproc)
            if iproc != Blockdata["main_node"]:
                dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
                orien_group_assignment[iproc_image_start:iproc_image_end] = dummy[:]
                del dummy
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    orien_group_assignment = sp_utilities.wrap_mpi_bcast(orien_group_assignment, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
    ptls_in_orien_groups = [None for iref in range(len(refa_vecs))]
    for iorien in range(len(refa_vecs)):
        if iorien % Blockdata["nproc"] == Blockdata["myid"]: ptls_in_orien_groups[iorien] = findall_dict(iorien,
                                                                                                         orien_group_assignment)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    for iorien in range(len(refa_vecs)):
        if iorien % Blockdata["nproc"] != Blockdata["main_node"]:
            if iorien % Blockdata["nproc"] == Blockdata["myid"]: sp_utilities.wrap_mpi_send(ptls_in_orien_groups[iorien],
                                                                               Blockdata["main_node"], mpi.MPI_COMM_WORLD)
            if Blockdata["myid"] == Blockdata["main_node"]: ptls_in_orien_groups[iorien] = sp_utilities.wrap_mpi_recv(
                iorien % Blockdata["nproc"], mpi.MPI_COMM_WORLD)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    zero_member_group_found = 0

    if Blockdata["myid"] == Blockdata["main_node"]:
        small_groups = []
        for iorien in range(len(refa_vecs)):
            if len(ptls_in_orien_groups[iorien]) < Tracker["min_orien_group_size"]: small_groups.append(iorien)
            if len(ptls_in_orien_groups[iorien]) == 0:
                zero_member_group_found += 1
        matched_pairs = find_neighborhood(refa_vecs, small_groups)
        if len(matched_pairs) >= 1:
            ptls_in_orien_groups = reassign_ptls_in_orien_groups(ptls_in_orien_groups, matched_pairs)
    else:
        ptls_in_orien_groups = 0
    zero_member_group_found = sp_utilities.bcast_number_to_all(zero_member_group_found, Blockdata["main_node"], mpi.MPI_COMM_WORLD)
    ptls_in_orien_groups = sp_utilities.wrap_mpi_bcast(ptls_in_orien_groups, Blockdata["main_node"], mpi.MPI_COMM_WORLD)

    del refa_vecs, refa
    del local_orien_group_assignment
    del data_angles
    del orien_group_assignment
    del sym_class
    return ptls_in_orien_groups














































































































































def copy_refinement_tracker(tracker_refinement):
    global Tracker, Blockdata
    for key, value in Tracker:
        try:
            value_refinement = tracker_refinement[key]
            if value == None and value_refinement != None: Tracker[key] = value_refinement
        except:
            if Blockdata["myid"] == Blockdata["main_node"]: print(key, " in sorting set as ", value,
                                                                  ", while in refinement, it is set as ",
                                                                  value_refinement)
    return


def print_dict(dict, theme):
    spaces = "                    "
    exclude = ["nodes", "yr", "output", "shared_comm", "bckgnoise", "myid", "myid_on_node", "accumulatepw", \
               "chunk_dict", "PW_dict", "full_list", "rshifts", "refang", "chunk_dict", "PW_dict", "full_list",
               "rshifts", \
               "refang", "sorting_data_list", "partition", "constants", "random_assignment"]
    for key, value in sorted(dict.items()):
        pt = True
        for ll in exclude:
            if (key == ll):
                pt = False
                break
        if pt:  print("                    => ", key + spaces[len(key):], ":  ", value)


#
# - "Tracker" (dictionary) object
#   Keeps the current state of option settings and dataset
#   (i.e. particle stack, reference volume, reconstructed volume, and etc)
#   Each iteration is allowed to add new fields/keys
#   if necessary. This happes especially when type of 3D Refinement or metamove changes.
#   Conceptually, each iteration will be associated to a specific Tracker state.
#   Therefore, the list of Tracker state represents the history of process.
#
#   This can be used to restart process from an arbitrary iteration.
#                     rec3d for sorting








































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































def do3d(procid, data, newparams, refang, rshifts, norm_per_particle, myid, mpi_comm=-1):
    global Tracker, Blockdata
    #  Without filtration
    pass#IMPORTIMPORTIMPORT from sp_reconstruction import recons3d_trl_struct_MPI
    if (mpi_comm < -1): mpi_comm = MPI_COMM_WORDLD
    if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
        if (procid == 0):
            if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")): os.mkdir(
                os.path.join(Tracker["directory"], "tempdir"))
    shrinkage = float(Tracker["nxinit"]) / float(Tracker["constants"]["nnxo"])
    tvol, tweight, trol = sp_reconstruction.recons3d_trl_struct_MPI(myid=Blockdata["subgroup_myid"], main_node=Blockdata["main_node"],
                                                  prjlist=data, \
                                                  paramstructure=newparams, refang=refang,
                                                  rshifts_shrank=[[q[0] * shrinkage, q[1] * shrinkage] for q in
                                                                  rshifts], \
                                                  delta=Tracker["delta"], CTF=Tracker["constants"]["CTF"],
                                                  upweighted=False, mpi_comm=mpi_comm, \
                                                  target_size=(2 * Tracker["nxinit"] + 3),
                                                  avgnorm=Tracker["avgvaradj"][procid],
                                                  norm_per_particle=norm_per_particle)
    if Blockdata["subgroup_myid"] == Blockdata["main_node"]:
        tvol.set_attr("is_complex", 0)
        tvol.write_image(
            os.path.join(Tracker["directory"], "tempdir", "tvol_%01d_%03d.hdf" % (procid, Tracker["mainiteration"])))
        tweight.write_image(
            os.path.join(Tracker["directory"], "tempdir", "tweight_%01d_%03d.hdf" % (procid, Tracker["mainiteration"])))
        trol.write_image(
            os.path.join(Tracker["directory"], "tempdir", "trol_%01d_%03d.hdf" % (procid, Tracker["mainiteration"])))
    mpi.mpi_barrier(mpi_comm)
    return


#######
def do3d_sorting_groups_rec3d(iteration, masterdir, log_main):
    global Tracker, Blockdata
    pass#IMPORTIMPORTIMPORT from sp_utilities import get_im
    # reconstruct final unfiltered volumes from sorted clusters
    keepgoing = 1
    ### ====
    fsc143 = 0
    fsc05 = 0
    Tracker["fsc143"] = 0
    Tracker["fsc05"] = 0
    res_05 = Tracker["number_of_groups"] * [0]
    res_143 = Tracker["number_of_groups"] * [0]
    Tracker["directory"] = masterdir
    Tracker["constants"]["masterdir"] = masterdir
    Tracker["maxfrad"] = Tracker["nxinit"] // 2
    ####
    if Blockdata["no_of_groups"] > 1:
        sub_main_node_list = [-1 for i in range(Blockdata["no_of_groups"])]
        for index_of_colors in range(Blockdata["no_of_groups"]):
            for iproc in range(Blockdata["nproc"] - 1):
                if Blockdata["myid"] == iproc:
                    if Blockdata["color"] == index_of_colors and Blockdata["myid_on_node"] == 0:
                        sub_main_node_list[index_of_colors] = Blockdata["myid"]
                    sp_utilities.wrap_mpi_send(sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
                if Blockdata["myid"] == Blockdata["last_node"]:
                    dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
                    for im in range(len(dummy)):
                        if dummy[im] > -1: sub_main_node_list[im] = dummy[im]
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        sp_utilities.wrap_mpi_bcast(sub_main_node_list, Blockdata["last_node"], mpi.MPI_COMM_WORLD)
        ####
        if Tracker["number_of_groups"] % Blockdata["no_of_groups"] == 0:
            nbig_loop = Tracker["number_of_groups"] // Blockdata["no_of_groups"]
        else:
            nbig_loop = Tracker["number_of_groups"] // Blockdata["no_of_groups"] + 1

        big_loop_colors = [[] for i in range(nbig_loop)]
        big_loop_groups = [[] for i in range(nbig_loop)]
        nc = 0
        while nc < Tracker["number_of_groups"]:
            im = nc // Blockdata["no_of_groups"]
            jm = nc % Blockdata["no_of_groups"]
            big_loop_colors[im].append(jm)
            big_loop_groups[im].append(nc)
            nc += 1
        #####
        for iloop in range(nbig_loop):
            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                Clusterdir = os.path.join(Tracker["directory"], "Cluster%d" % index_of_group, "main%03d" % iteration)
                if (Blockdata["myid"] == Blockdata["last_node"]):
                    tvol2 = sp_utilities.get_im(os.path.join(Clusterdir, "tempdir", "tvol_0_%03d.hdf" % iteration))
                    tweight2 = sp_utilities.get_im(os.path.join(Clusterdir, "tempdir", "tweight_0_%03d.hdf" % iteration))
                    treg2 = sp_utilities.get_im(os.path.join(Clusterdir, "tempdir", "trol_0_%03d.hdf" % iteration))
                    tag = 7007
                    sp_utilities.send_EMData(tvol2, sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
                    sp_utilities.send_EMData(tweight2, sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
                    sp_utilities.send_EMData(treg2, sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
                elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
                    tag = 7007
                    tvol2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
                    tweight2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
                    treg2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                if Blockdata["color"] == index_of_colors:
                    if (Blockdata["myid_on_node"] != 0):
                        tvol2 = sp_utilities.model_blank(1)
                        tweight2 = sp_utilities.model_blank(1)
                        treg2 = sp_utilities.model_blank(1)
                    tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False,
                                        color=index_of_colors)  # has to be False!!!
                    del tweight2, treg2
                mpi.mpi_barrier(Blockdata["shared_comm"])
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
                    tag = 7007
                    sp_utilities.send_EMData(tvol2, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
                elif (Blockdata["myid"] == Blockdata["last_node"]):
                    tag = 7007
                    tvol2 = sp_utilities.recv_EMData(sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
                    tvol2.write_image(
                        os.path.join(Tracker["directory"], "vol_unfiltered_0_grp%03d.hdf" % index_of_group))
                    del tvol2
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                Clusterdir = os.path.join(Tracker["directory"], "Cluster%d" % index_of_group, "main%03d" % iteration)
                if (Blockdata["myid"] == Blockdata["last_node"]):
                    tvol2 = sp_utilities.get_im(os.path.join(Clusterdir, "tempdir", "tvol_1_%03d.hdf" % iteration))
                    tweight2 = sp_utilities.get_im(os.path.join(Clusterdir, "tempdir", "tweight_1_%03d.hdf" % iteration))
                    treg2 = sp_utilities.get_im(os.path.join(Clusterdir, "tempdir", "trol_1_%03d.hdf" % iteration))
                    tag = 7007
                    sp_utilities.send_EMData(tvol2, sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
                    sp_utilities.send_EMData(tweight2, sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
                    sp_utilities.send_EMData(treg2, sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)

                elif (Blockdata["myid"] == sub_main_node_list[index_of_colors]):
                    tag = 7007
                    tvol2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
                    tweight2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
                    treg2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                Tracker["maxfrad"] = Tracker["nxinit"] // 2
                if Blockdata["color"] == index_of_colors:
                    if (Blockdata["myid_on_node"] != 0):
                        tvol2 = sp_utilities.model_blank(1)
                        tweight2 = sp_utilities.model_blank(1)
                        treg2 = sp_utilities.model_blank(1)
                    tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False,
                                        color=index_of_colors)  # has to be False!!!
                    del tweight2, treg2
                mpi.mpi_barrier(Blockdata["shared_comm"])
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            for im in range(len(big_loop_colors[iloop])):
                index_of_group = big_loop_groups[iloop][im]
                index_of_colors = big_loop_colors[iloop][im]
                if (Blockdata["color"] == index_of_colors) and (Blockdata["myid_on_node"] == 0):
                    tag = 7007
                    sp_utilities.send_EMData(tvol2, Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
                elif (Blockdata["myid"] == Blockdata["last_node"]):
                    tag = 7007
                    tvol2 = sp_utilities.recv_EMData(sub_main_node_list[index_of_colors], tag, mpi.MPI_COMM_WORLD)
                    tvol2.write_image(
                        os.path.join(Tracker["directory"], "vol_unfiltered_1_grp%03d.hdf" % index_of_group))
                    del tvol2
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    else:
        Tracker["maxfrad"] = Tracker["nxinit"] // 2
        for index_of_group in range(Tracker["number_of_groups"]):
            Clusterdir = os.path.join(Tracker["directory"], "Cluster%d" % index_of_group, "main%03d" % iteration)

            if (Blockdata["myid"] == Blockdata["last_node"]):
                tvol2 = sp_utilities.get_im(os.path.join(Clusterdir, "tempdir", "tvol_0_%03d.hdf" % iteration))
                tweight2 = sp_utilities.get_im(os.path.join(Clusterdir, "tempdir", "tweight_0_%03d.hdf" % iteration))
                treg2 = sp_utilities.get_im(os.path.join(Clusterdir, "tempdir", "trol_0_%03d.hdf" % iteration))
                tag = 7007
                sp_utilities.send_EMData(tvol2, Blockdata["main_node"], tag, mpi.MPI_COMM_WORLD)
                sp_utilities.send_EMData(tweight2, Blockdata["main_node"], tag, mpi.MPI_COMM_WORLD)
                sp_utilities.send_EMData(treg2, Blockdata["main_node"], tag, mpi.MPI_COMM_WORLD)
            elif (Blockdata["myid"] == Blockdata["main_node"]):
                tag = 7007
                tvol2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
                tweight2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
                treg2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            if Blockdata["myid"] != Blockdata["main_node"]:
                tvol2 = sp_utilities.model_blank(1)
                tweight2 = sp_utilities.model_blank(1)
                treg2 = sp_utilities.model_blank(1)
            tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False, color=0)  # has to be False!!!
            del tweight2, treg2
            if (Blockdata["myid"] == Blockdata["main_node"]):
                tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_0_grp%03d.hdf" % index_of_group))
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            if (Blockdata["myid"] == Blockdata["last_node"]):
                tvol2 = sp_utilities.get_im(os.path.join(Clusterdir, "tempdir", "tvol_1_%03d.hdf" % iteration))
                tweight2 = sp_utilities.get_im(os.path.join(Clusterdir, "tempdir", "tweight_1_%03d.hdf" % iteration))
                treg2 = sp_utilities.get_im(os.path.join(Clusterdir, "tempdir", "trol_1_%03d.hdf" % iteration))
                tag = 7007
                sp_utilities.send_EMData(tvol2, Blockdata["main_node"], tag, mpi.MPI_COMM_WORLD)
                sp_utilities.send_EMData(tweight2, Blockdata["main_node"], tag, mpi.MPI_COMM_WORLD)
                sp_utilities.send_EMData(treg2, Blockdata["main_node"], tag, mpi.MPI_COMM_WORLD)
            elif (Blockdata["myid"] == Blockdata["main_node"]):
                tag = 7007
                tvol2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
                tweight2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
                treg2 = sp_utilities.recv_EMData(Blockdata["last_node"], tag, mpi.MPI_COMM_WORLD)
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            if Blockdata["myid"] != Blockdata["main_node"]:
                tvol2 = sp_utilities.model_blank(1)
                tweight2 = sp_utilities.model_blank(1)
                treg2 = sp_utilities.model_blank(1)
            tvol2 = steptwo_mpi(tvol2, tweight2, treg2, None, False, color=0)  # has to be False!!!
            del tweight2, treg2
            if (Blockdata["myid"] == Blockdata["main_node"]):
                tvol2.write_image(os.path.join(Tracker["directory"], "vol_unfiltered_1_grp%03d.hdf" % index_of_group))
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    keepgoing = sp_utilities.bcast_number_to_all(keepgoing, source_node=Blockdata["main_node"],
                                    mpi_comm=mpi.MPI_COMM_WORLD)  # always check
    Tracker = sp_utilities.wrap_mpi_bcast(Tracker, Blockdata["main_node"])
    if keepgoing == 0: sp_global_def.ERROR("do3d_sorting_groups_trl_iter  %s" % os.path.join(Tracker["directory"], "tempdir"),
                             "do3d_sorting_groups_trl_iter", 1, Blockdata["myid"])
    return


####=====<--------
### nofsc rec3d























































































































































































































































































































































































































































































































































































































































































































































































































def do3d_sorting_group_insertion_smearing(sdata, sparamstructure, snorm_per_particle, randomset=2):
    global Tracker, Blockdata
    if (Blockdata["myid"] == Blockdata["nodes"][0]):
        if not os.path.exists(os.path.join(Tracker["directory"], "tempdir")): os.mkdir(
            os.path.join(Tracker["directory"], "tempdir"))
    if randomset == 1:
        for index_of_groups in range(Tracker["number_of_groups"]):
            for procid in range(2, 3):
                tvol, tweight, trol = recons3d_trl_struct_group_MPI(myid=Blockdata["myid"],
                                                                    main_node=Blockdata["nodes"][procid // 2],
                                                                    prjlist=sdata, random_subset=procid,
                                                                    group_ID=index_of_groups,
                                                                    paramstructure=sparamstructure,
                                                                    norm_per_particle=snorm_per_particle,
                                                                    CTF=Tracker["constants"]["CTF"],
                                                                    mpi_comm=None, upweighted=False,
                                                                    target_size=(2 * Tracker["nxinit"] + 3),
                                                                    nosmearing=Tracker["nosmearing"])
                if (Blockdata["myid"] == Blockdata["nodes"][procid // 2]):
                    tvol.set_attr("is_complex", 0)
                    tvol.write_image(
                        os.path.join(Tracker["directory"], "tempdir", "tvol_%d_%d.hdf" % (procid, index_of_groups)))
                    tweight.write_image(
                        os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf" % (procid, index_of_groups)))
                    trol.write_image(
                        os.path.join(Tracker["directory"], "tempdir", "trol_%d_%d.hdf" % (procid, index_of_groups)))
                del tvol
                del tweight
                del trol
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    else:
        for index_of_groups in range(Tracker["number_of_groups"]):
            for procid in range(3):
                tvol, tweight, trol = recons3d_trl_struct_group_MPI(myid=Blockdata["myid"],
                                                                    main_node=Blockdata["nodes"][procid // 2],
                                                                    prjlist=sdata, random_subset=procid,
                                                                    group_ID=index_of_groups,
                                                                    paramstructure=sparamstructure,
                                                                    norm_per_particle=snorm_per_particle,
                                                                    CTF=Tracker["constants"]["CTF"],
                                                                    mpi_comm=None, upweighted=False,
                                                                    target_size=(2 * Tracker["nxinit"] + 3),
                                                                    nosmearing=Tracker["nosmearing"])
                if (Blockdata["myid"] == Blockdata["nodes"][procid // 2]):
                    tvol.set_attr("is_complex", 0)
                    tvol.write_image(
                        os.path.join(Tracker["directory"], "tempdir", "tvol_%d_%d.hdf" % (procid, index_of_groups)))
                    tweight.write_image(
                        os.path.join(Tracker["directory"], "tempdir", "tweight_%d_%d.hdf" % (procid, index_of_groups)))
                    trol.write_image(
                        os.path.join(Tracker["directory"], "tempdir", "trol_%d_%d.hdf" % (procid, index_of_groups)))
                del tvol
                del tweight
                del trol
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    return


### rec3d
####=====<-------MEM related functions
























def resident(since=0.0):
    '''Return resident memory usage in bytes.
    '''
    return _VmB('VmRSS:') - since


def stacksize(since=0.0):
    '''Return stack size in bytes.
    '''
    return _VmB('VmStk:') - since


#####==========-------------------------Functions for post processing




















































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































