import math
import numpy
# def INTERPOL_WRAP(i , n):
#     if i >= 0 :
#         y = (i) % (n)
#         return y
#     else:
#         if i< 0:
#             try:
#                 y = (n) - (-(i) % (n))
#             except:
#                 y = -(i) % (n)
#             return y
#         return y



def INTERPOL_WRAP(i, n):
    value = 0
    if i >= 0 :
        value = (i) % (n)
    elif -(i) % (n) :
        value = (n) - (-(i) % (n))
    else:
        value = 0
    return  value


def cubicXY(img, x , y):
    xi = int(numpy.floor(x))
    yi = int(numpy.floor(y))
    xi_n1 = xi - 1
    yi_n1 = yi - 1
    xi_p1 = xi + 1
    yi_p1 = yi + 1
    xi_p2 = xi + 2
    yi_p2 = yi + 2
    xf = x - xi
    yf = y - yi

    img_xdim = img.shape[0]
    img_ydim = img.shape[1]
    xi_n1 = INTERPOL_WRAP(xi_n1, img_xdim)
    yi_n1 = INTERPOL_WRAP(yi_n1, img_ydim)
    xi = INTERPOL_WRAP(xi, img_xdim)
    yi = INTERPOL_WRAP(yi, img_ydim)
    xi_p1 = INTERPOL_WRAP(xi_p1, img_xdim)
    yi_p1 = INTERPOL_WRAP(yi_p1, img_ydim)
    xi_p2 = INTERPOL_WRAP(xi_p2, img_xdim)
    yi_p2 = INTERPOL_WRAP(yi_p2, img_ydim)

    f00 = img[yi_n1][xi_n1]
    f01 = img[yi_n1][xi]
    f02 = img[yi_n1][xi_p1]
    f03 = img[yi_n1][xi_p2]

    f10 = img[yi][xi_n1]
    f11 = img[yi][xi]
    f12 = img[yi][xi_p1]
    f13 = img[yi][xi_p2]

    f20 = img[yi_p1][xi_n1]
    f21 = img[yi_p1][xi]
    f22 = img[yi_p1][xi_p1]
    f23 = img[yi_p1][xi_p2]

    f30 = img[yi_p2][xi_n1]
    f31 = img[yi_p2][xi]
    f32 = img[yi_p2][xi_p1]
    f33 = img[yi_p2][xi_p2]

    A = numpy.array([ -1.0/2.0,  3.0/2.0, -3.0/2.0,  1.0/2.0,
             1.0,     -5.0/2.0,  2.0,     -1.0/2.0,
             -1.0/2.0,  0.0,      1.0/2.0,  0.0,
            0.0,      1.0,      0.0,      0.0]).reshape(4,4)

    V = numpy.array([f00, f10, f20, f30,
                  f01, f11, f21, f31,
                  f02, f12, f22, f32,
                   f03, f13, f23, f33]).reshape(4,4)

    # A = numpy.transpose(A)
    # V = numpy.transpose(V)

    At = A
    At = numpy.transpose(At)

    AVA = numpy.matmul(numpy.matmul(A , V) , At)
    xx = numpy.array([xf * xf * xf, xf * xf, xf, 1.0])
    yy = numpy.array([yf * yf * yf, yf * yf, yf, 1.0])

    return xx.dot(numpy.inner(AVA, yy ))




def cubicXYgrad(img, x , y):
    xi = int(numpy.floor(x))
    yi = int(numpy.floor(y))
    xi_n1 = xi - 1
    yi_n1 = yi - 1
    xi_p1 = xi + 1
    yi_p1 = yi + 1
    xi_p2 = xi + 2
    yi_p2 = yi + 2
    xf = x - xi
    yf = y - yi

    img_xdim = img.shape[0]
    img_ydim = img.shape[1]
    xi_n1 = INTERPOL_WRAP(xi_n1, img_xdim)
    yi_n1 = INTERPOL_WRAP(yi_n1, img_ydim)
    xi = INTERPOL_WRAP(xi, img_xdim)
    yi = INTERPOL_WRAP(yi, img_ydim)
    xi_p1 = INTERPOL_WRAP(xi_p1, img_xdim)
    yi_p1 = INTERPOL_WRAP(yi_p1, img_ydim)
    xi_p2 = INTERPOL_WRAP(xi_p2, img_xdim)
    yi_p2 = INTERPOL_WRAP(yi_p2, img_ydim)

    f00 = img[yi_n1][xi_n1]
    f01 = img[yi_n1][xi]
    f02 = img[yi_n1][xi_p1]
    f03 = img[yi_n1][xi_p2]

    f10 = img[yi][xi_n1]
    f11 = img[yi][xi]
    f12 = img[yi][xi_p1]
    f13 = img[yi][xi_p2]

    f20 = img[yi_p1][xi_n1]
    f21 = img[yi_p1][xi]
    f22 = img[yi_p1][xi_p1]
    f23 = img[yi_p1][xi_p2]

    f30 = img[yi_p2][xi_n1]
    f31 = img[yi_p2][xi]
    f32 = img[yi_p2][xi_p1]
    f33 = img[yi_p2][xi_p2]

    A = numpy.array([ -1.0/2.0,  3.0/2.0, -3.0/2.0,  1.0/2.0,
             1.0,     -5.0/2.0,  2.0,     -1.0/2.0,
             -1.0/2.0,  0.0,      1.0/2.0,  0.0,
            0.0,      1.0,      0.0,      0.0]).reshape(4,4)

    V = numpy.array([f00, f10, f20, f30,
                  f01, f11, f21, f31,
                  f02, f12, f22, f32,
                   f03, f13, f23, f33]).reshape(4,4)

    # A = numpy.transpose(A)
    # V = numpy.transpose(V)

    At = A
    At = numpy.transpose(At)

    AVA = numpy.matmul(numpy.matmul(A , V) , At)
    xx = numpy.array([xf * xf * xf, xf * xf, xf, 1.0])
    yy = numpy.array([yf * yf * yf, yf * yf, yf, 1.0])

    xxd = numpy.array([3.0*xf*xf, 2.0*xf, 1.0, 0.0])
    yyd = numpy.array([3.0*yf*yf, 2.0*yf, 1.0, 0.0])

    return xxd.dot(numpy.inner(AVA, yy )) ,  xx.dot(numpy.inner(AVA, yyd ))















# img = numpy.arange(129600).reshape(360,360)
#
# epf =  cubicXY( img, -0.963542999999 , 0.81565000000)
#
#
# print('heloo')