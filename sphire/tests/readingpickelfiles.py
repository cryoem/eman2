from __future__ import print_function
import numpy
import pickle
import os
import EMAN2_cppwrap as e2cpp

ABSOLUTE_PATH =  os.path.dirname(os.path.realpath(__file__))
filepath = os.path.join(ABSOLUTE_PATH, "files/PICKLE.ornq")

# with open(filepath, 'rb') as rb:
#       (image,crefim,xrng,yrng,step,mode,numr,cnx,cny,deltapsi) = pickle.load(rb)
#       print(numpy.shape(image.get_3dview()))
#       print(numpy.shape(crefim.get_3dview()))
#       print(xrng)
#       print(yrng)
#       print(step)
#       print(mode)
#       print(numpy.shape(numr))
#       print(cnx)
#       print(cny)
#       print(deltapsi)

import io
import cPickle
import sys
import types
from functools import wraps
import fcntl

funclist = ['ali2d_single_iter', 'ringwe' ,'ormq_fast', 'prepref', 'prepare_refrings', 'proj_ali_incore',\
            'proj_ali_incore_local', 'ali_vol_func', 'align2d', 'align2d_scf', 'multalign2d_scf', \
            'parabl', 'shc', 'search_range', 'generate_list_of_reference_angles_for_search']

# #
# FUNCLIST = ['ali2d_MPI', 'ali2d_base', 'mref_ali3d_MPI', 'Kmref_ali3d_MPI', 'cpy', 'project3d', 'ali_vol',\
#             'recons3d_n_trl_MPI_one_node', 'pca', 'prepare_2d_forPCA', 'header', 'refvol', \
#              'within_group_refinement','ali3d_mref_Kmeans_MPI', 'mref_ali3d_EQ_Kmeans']

FUNCLIST = ['ormq_fast','ringwe' , 'prepref','proj_ali_incore','proj_ali_incore_local']


#sparx_multi_shc
FUNCTLIST = ['orient_params', 'find_common_subset', 'ali3d_multishc', 'ali3d_multishc_2', 'multi_shc',\
            'mirror_and_reduce_dsym', 'do_volume' ]

#sparx_pixe_error
FUNCLIST = ['pixel_error_2D','max_3D_pixel_error','angle_ave','angle_diff','angle_diff_sym', 'align_diff_params',\
            'multi_align_stability']


#sparx_projections
FUNCLIST = ['project', 'prgs', 'prgl', 'prgq', 'prg', 'prep_vol', 'gen_rings_ctf']


#sparx_reconstruction
FUNCLIST = ['insert_slices', 'insert_slices_pdf', 'recons3d_4nn_MPI', 'recons3d_trl_struct_MPI', 'recons3d_4nn_ctf', \
            'recons3d_4nn_ctf_MPI', 'recons3d_nn_SSNR_MPI', 'prepare_recons', 'prepare_recons_ctf', 'recons_from_fftvol', \
            'recons_ctf_from_fftvol', 'get_image_size', 'rec3D_MPI', 'rec3D_MPI_noCTF', 'prepare_recons_ctf_two_chunks',\
            'rec3D_two_chunks_MPI']

#sparx_statistics  (functions inside the classes are not included)
FUNCLIST = [ 'add_ave_varf_MPI', 'sum_oe', 'ave_var', 'ave_series', 'varf2d_MPI', 'varf3d_MPI', 'ccc', \
            'fsc', 'fsc_mask', 'locres', 'histogram', 'k_means_match_clusters_asg_new', 'hist_list', \
            'linreg', 'pearson', 'table_stat', 'mono', 'k_means_stab_bbenum', 'k_means_match_bbenum', 'scale_fsc_datasetsize' ]


#sparx_user_functions   (functions inside the classes are not included)
FUNCLIST = ['ref_ali2d', 'ref_ali2d_c', 'julien', 'ref_ali2d_m', 'ref_ali3dm', 'ref_sort3d', 'ref_ali3dm_ali_50S',\
            'ref_random', 'ref_ali3d', 'helical', 'helical2', 'reference3', 'reference4', 'ref_aliB_cone', \
            'ref_7grp', 'spruce_up', 'spruce_up_variance', 'minfilt', 'ref_ali3dm_new', 'spruce_up_var_m', \
            'steady', 'constant', 'temp_dovolume', 'dovolume', 'do_volume_mask', 'do_volume_mrk02', 'do_volume_mrk03', \
            'do_volume_mrk04', 'do_volume_mrk05', 'build_user_function']



#sparx_utilities (function within a function is not included, functions inside the classes are not included)
FUNCLIST = [ 'amoeba', 'compose_transform2', 'compose_transform3', 'combine_params2', 'inverse_transform2', \
             'drop_image', 'drop_spider_doc', 'even_angles', 'even_angles_cd', 'find', 'gauss_edge', \
             'get_image', 'get_im', 'get_image_data', 'get_symt', 'get_input_from_string', 'model_circle', \
            'model_gauss', 'model_gauss_noise', 'model_blank', 'peak_search', 'print_list_format', 'pad', \
             'chooseformat', 'read_text_row', 'write_text_row', 'read_text_file', 'write_text_file',  'rotate_shift_params', \
             'reshape_1d', 'estimate_3D_center_MPI', 'rotate_3D_shift', 'set_arb_params', 'get_arb_params', 'reduce_EMData_to_root', \
             'bcast_compacted_EMData_all_to_all', 'gather_compacted_EMData_to_root', 'bcast_EMData_to_all', 'send_EMData', \
             'recv_EMData', 'bcast_number_to_all', 'bcast_list_to_all', 'recv_attr_dict', 'send_attr_dict', \
             'recv_attr_dict_bdb', 'print_begin_msg', 'print_end_msg', 'print_msg', 'read_fsc', 'circumference', \
             'write_headers', 'write_header', 'file_type', 'get_params2D', 'set_params2D', 'get_params3D', 'set_params3D', \
             'get_params_proj', 'set_params_proj', 'get_ctf', 'same_ctf', 'generate_ctf', 'delete_bdb', \
             'disable_bdb_cache', 'getvec', 'getfvec', 'nearest_fang', 'nearest_many_full_k_projangles', 'assign_projdirs_f', \
             'angles_to_normals', 'angular_occupancy', 'angular_histogram', 'balance_angular_distribution', 'symmetry_neighbors', \
             'rotation_between_anglesets', 'angle_between_projections_directions', 'get_pixel_size', 'set_pixel_size', \
             'lacos', 'nearest_proj', 'findall', 'pack_message', 'unpack_message', 'update_tag', 'wrap_mpi_send', 'wrap_mpi_recv', \
             'wrap_mpi_bcast', 'wrap_mpi_gatherv', 'get_colors_and_subsets', 'wrap_mpi_split', 'get_dist', 'eliminate_moons', \
             'combinations_of_n_taken_by_k', 'cmdexecute', 'string_found_in_file', 'get_latest_directory_increment_value', \
             'if_error_then_all_processes_exit_program', 'get_shrink_data_huang', 'getindexdata',  \
             'store_value_of_simple_vars_in_json_file', 'convert_json_fromunicode', 'get_sorting_attr_stack', \
             'get_sorting_params_refine', 'parsing_sorting_params', 'fill_in_mpi_list', 'sample_down_1D_curve', \
             'get_initial_ID', 'print_upper_triangular_matrix', 'convertasi', 'prepare_ptp', 'print_dict',  \
             'get_resolution_mrk01', 'partition_to_groups', 'partition_independent_runs', 'merge_groups', 'save_alist', \
             'margin_of_error', 'do_two_way_comparison', 'select_two_runs', 'counting_projections', 'unload_dict', \
             'load_dict', 'get_stat_proj', 'create_random_list', 'recons_mref', 'apply_low_pass_filter', 'get_groups_from_partition', \
             'get_complementary_elements', 'update_full_dict', 'count_chunk_members', 'remove_small_groups', 'get_number_of_groups', \
             'angular_distribution', 'tabessel']


for entry in FUNCLIST[:]:
      print(entry)

# def pickle_arguments(f):
#     # @wraps(f)
#     def decorated(*args,**kwargs):
#         global FUNCLIST
#         for entry in FUNCLIST[:]:
#             if f.__name__ == entry:
#                 if os.path.isfile('alignment'+'.'+ f.__name__):
#                     pass
#                     print('file already exists')
#                 else:
#                     with open('alignment'+'.'+ f.__name__, 'wb') as wb:
#                         cPickle.dump((args,kwargs),wb)
#                     wb.close()
#                     print('alignment'+'.'+ f.__name__)
#                     print('success')
#                     FUNCLIST.remove(entry)
#         print('End of the section:')
#         return f(*args, **kwargs)
#     return decorated




def pickle_arguments(f):
    def decorated(*args,**kwargs):
        global FUNCLIST
        print('FUNCTION', f.__name__, 'START', file=sys.stderr)
        for entry in FUNCLIST[:]:
            if f.__name__ == entry:
                if os.path.isfile('alignment'+'.'+ f.__name__):
                    print("file already created",file=sys.stderr)
                    pass
                else:
                    print("writing into file",file=sys.stderr)
                    with open('alignment'+'.'+ f.__name__, 'wb') as wb:
                        cPickle.dump((args,kwargs),wb)
                        cPickle.dumps(Foo)
                    wb.close()
                    print('alignment'+'.'+ f.__name__,file=sys.stderr)
                    print('success')
                    FUNCLIST.remove(entry)
        print('FUNCTION', f.__name__, 'DONE', file=sys.stderr)
        return f(*args, **kwargs)
    return decorated



def prepref(a,b,c = 'string', dd = 0.2 , e= ''):
      c = a+b
      return(a, b , c )


@pickle_arguments
def ormq_fast(a,b,c = 'string', dd = 0.2):
      c = a+b
      return(a, b , c )


@pickle_arguments
def ringwe(a,b,c = 'string', dd = 0.2):
      c = a+b
      return(a, b , c )

a = prepref(2,6, 'dumm', 0.5, '')
print(a)

a = ormq_fast(3,4, 'dumm', 0.5)
print(a)

a = ringwe(7,1, 'dumm', 0.5)
print(a)

a = prepref(8,5, 'dumm', 0.5, '')
print(a)
