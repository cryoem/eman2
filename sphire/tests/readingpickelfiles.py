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



FUNCTLIST = ['orient_params', 'find_common_subset', 'ali3d_multishc', 'ali3d_multishc_2', 'multi_shc',\
            'mirror_and_reduce_dsym', 'do_volume' ]

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
class Foo:
    attr = 'A class attribute'
    Protected attributes:
        float alpha
        float bet



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
