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


import pickle
from functools import wraps

funclist = ['ali2d_single_iter', 'ringwe' ,'ormq_fast', 'prepref', 'prepare_refrings', 'proj_ali_incore',\
            'proj_ali_incore_local', 'ali_vol_func', 'align2d', 'align2d_scf', 'multalign2d_scf', \
            'parabl', 'shc', 'search_range', 'generate_list_of_reference_angles_for_search']

#
# funclist = ['ali2d_MPI','ali2d_base', 'mref_ali3d_MPI', 'Kmref_ali3d_MPI', 'cpy', 'project3d',\
#             'ali_vol', 'recons3d_n_trl_MPI_one_node', 'pca', 'prepare_2d_forPCA', 'header']

# def pickle_arguments(f):
#       @wraps(f)
#       def decorated(*args,**kwargs):
#             # print('alignment'+'.'+ f.__name__)
#             with open('alignment'+'.'+ f.__name__, 'wb') as wb:
#                   pickle.dump((args,kwargs),wb)
#             return f(*args,**kwargs)
#       return decorated
#
#
#
# def findfunc_in_list(f):
#       @wraps(f)
#       def decore(*args, **kwargs):
#             # print('Finding function if it is already written')
#             for entry in funclist:
#                   if f.__name__ == entry:
#                         # print("True statement executed")
#                         # print(entry)
#                         funclist[funclist.index(entry)] = ''
#                         return f(*args, **kwargs)
#                   else:
#                         # print("False Statement executed")
#                         pass
#             # return f(*args, **kwargs)
#       return decore

rmindex = -1
def pickle_arguments(f):
      @wraps(f)
      def decorated(*args,**kwargs):
            for entry in funclist:
                  global rmindex
                  if (f.__name__ == entry) and (funclist[funclist.index(entry)] != ''):
                        rmindex = funclist.index(entry)
                        with open('alignment'+'.'+ f.__name__, 'wb') as wb:
                              pickle.dump((args,kwargs),wb)
                        funclist[funclist.index(entry)] = ''
                        return f(*args,**kwargs)
                  elif (funclist[rmindex] == ''):
                        return f(*args, **kwargs)
      return decorated


# def findfunc_in_list(f):
#       @wraps(f)
#       def decore(*args, **kwargs):
#             for entry in funclist:
#                   if f.__name__ == entry:
#                         funclist[funclist.index(entry)] = ''
#                         return f(*args, **kwargs)
#                   else:
#                         pass
#       return decore

# @findfunc_in_list
@pickle_arguments
def prepref(a,b,c = 'string', dd = 0.2):
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

a = prepref(2,6, 'dumm', 0.5)
print(a)

a = ormq_fast(3,4, 'dumm', 0.5)
print(a)

a = ringwe(7,1, 'dumm', 0.5)
print(a)

a = prepref(8,5, 'dumm', 0.5)
print(a)