
######################################
#  The even case, 3D:  Part 1
######################################

N=4;  # This is the size of the original image
ti = EMData(N,N,N);
ti.to_zero();

tiB    = ti.process('xform.phaseorigin.tocorner')
tiBFFT = tiB.do_fft();

#  The resulting FFT should have
# nx=6, ny=4, nz=4

nx=tiBFFT.get_xsize();
ny=tiBFFT.get_ysize();
nz=tiBFFT.get_zsize();
print nx,ny,nz

# the 3 values of kx are 0, 1, 2
# the 4 values of ky  are 0, 1, 2, -1(3)
#   the above convention is SteveL's
#  get_complex_at will return 0, when ky = 3
# nx=6, because there are real and complex 
# values in this "direction"



#  So ny=nz=N, the side of the real image
#  And nx=N+2

# Let's set the values of the array

for kz in range(ny):
 for ky in range(ny):
  for kx in range(nx/2):
    tiBFFT.set_value_at(2*kx  ,ky, kz,  kx+ky*nx/2 +kz*ny*nx/2);
    tiBFFT.set_value_at(2*kx+1,ky, kz, -(kx+ky*nx/2+kz*ny*nx/2));
    
# Now the real and complex values are negatives
#     of one another
# the minimum magnitude is zero
# the maximum magnitude is 23

#  This does not correspond to a
#  real image, but that is no matter for this section



# Let's retrieve the information using get_complex_at
# kx =0, ky=0, kx+4 ky=0
tiBFFT.get_complex_at(0,0,0)  # gives 0 correctly

# kx =1, ky=1,, kx+4 ky=5
tiBFFT.get_complex_at(1,1);   #  gives 5 -5j correctly
tiBFFT.get_value_at(2,1); tiBFFT.get_value_at(3,1);  # is 5, -5

# kx =3, ky=-1 ( or 5),, kx+4 ky=23
tiBFFT.get_complex_at(3,5-6)  # should be 23 -23j
tiBFFT.get_value_at(6,5); tiBFFT.get_value_at(7,5);  # is 23, -23


#  Friedel pairs
# kx =-1, ky=-1,
#  This is Friedel pair of 
# kx=1 , ky=1,   kx+4*ky=5
tiBFFT.get_complex_at(-1,-1);   #  gives 5 + 5j correctly  => 5, -5
tiBFFT.get_value_at(2,1); tiBFFT.get_value_at(3,1);  # is 5, -5

tiBFFT.get_complex_at(-3,-1);   #  gives 5 + 5j correctly  => 5, -5
tiBFFT.get_value_at(2,1); tiBFFT.get_value_at(3,1);  # is 5, -5


###########################################
## Part 2.   Counting numbers of variables
#     get_complex_at   for N even, 3D

# Naively there are a total of 96=(N+2)NN real values in the complex array
#  that stores the 6 by 4 by 4 FFT image
# However the original image housed only   4*4*4 =64 real values

# The redundancies are 
# kx=(0 or 2), ky= (0 or 2), kz= (0 or 2)   value is  real        8


# [ (kx,ky) equal 0 or 2 ], kz=1,3  Friedel pairs     4 *2     
# [ (kx,kz) equal 0 or 2 ], ky=1,3  Friedel pairs     4 *2
# [ kx =0 or 2, and (ky,kz) = (1,1) or (3,3)]  Friedel pairs     2*2
# [ kx =0 or 2, and (ky,kz) = (1,3) or (3,1)]  Friedel pairs     2*2

# Generally there are  8 reals
#  Friedel pair redundancies 8 *(N-2) of the type of the first two lines of last paragraph
#  Friedel pairs redundancies   2*(N-2)*(N-2)   of the last two lines of last paragraph

# so there are total of  32 conditions = 2 N^2
#  Let's make sure they hold


N=6
ti=EMData(N,N); ti.to_zero();
for jx in range(N):
  for jy in range(N):
    ti.set_value_at(jx,jy, random());

tiB    = ti.process('xform.phaseorigin.tocorner')
tiBFFT = tiB.do_fft();


#  4 values
print tiBFFT.get_complex_at(0,0);     #  gives real value
print tiBFFT.get_complex_at(0,N/2);   #  gives real value
print tiBFFT.get_complex_at(N/2,0);   #  gives real value
print tiBFFT.get_complex_at(N/2,N/2); #  gives real value


#  2(N/2 -1) 2 redundancies
for kx in [0, N/2]:#    so 2 values
  for ky in range(1,N/2): #  (N/2-1) values
    print tiBFFT.get_complex_at(kx,ky);
    print tiBFFT.get_complex_at(kx,-ky);# this is the cc of above line
    print 'Last two values should be complex conjugates';

# so real values = complex values
# N^2= (N+2)(N) -4   -(2)(N/2-1)(2) 

###########################################################
#   Part 3  set_complex_at
#
#  So there are  2N special points to worry about when using set_complex_at
#  4  of them we must make sure are real
#  2N-4 of them we must set the Friedel pair

#   This is tricky, because  the set command may influence
#    multiple Points
#  So if we set anything with kx=0,3 or ky=0,3, (that is 0 or N/2)
#      we have to ensure that something else is true
#   in order that our complex image is an FFT of a real image
#  (the Friedel pair must be a complex conjugate)

# set_complex_at will accept those parts of the set that are admissible
#  For examples

tiBFFT.set_complex_at(3,0,1+2j)
tiBFFT.get_complex_at(3,0)#  gives 1

# set_complex_at won't do wrapping 
tiBFFT.set_complex_at(3,4,1+2j)

tiBFFT.get_complex_at(3,4)# Out[62]: 0j
tiBFFT.get_complex_at(3,-2) # (1.64420747756958-0.003193262964487076j)

# the redundant information it will create correctly

tiBFFT.set_complex_at(0,1,1+2j)
tiBFFT.get_complex_at(0,1)#   1 +2j
tiBFFT.get_complex_at(0,-1)#  1-2j


# Let's check this completely
for kx in [0, N/2]:#    so 2 values
  for ky in range(1,N/2):#  (N/2-1) values
    tiBFFT.set_complex_at(kx,ky,kx);
    print kx,tiBFFT.get_complex_at(kx,-ky);# this is the cc of above line
    print 'Last two values should be complex conjugates';


#############################################

######################################
#  The odd case, 2D:  Part 1
######################################

N=5;  # This is the size of the original image
ti = EMData(N,N);
ti.to_zero();

tiB    = ti.process('xform.phaseorigin.tocorner')
tiBFFT = tiB.do_fft();

#  The resulting FFT should have
# nx=6, ny=5

nx=tiBFFT.get_xsize();
ny=tiBFFT.get_ysize();
print nx,ny

# the 3 values of kx are 0, 1, 2
# the 5 values of ky  are 0, 1, 2,  -2(3), -1(4)
#   the above convention is SteveL's
#  get_complex_at will return 0, when ky = 3
# nx=6, because there are real and complex 
# values in this "direction"



#  So ny=N, the side of the real image
#  And nx=N+1 (as opposed to nx=N+2)
# General formula is nx=  2*ceil((N+1)/2)

# Let's set the values of the array

for ky in range(ny):
  for kx in range(nx/2):
    tiBFFT.set_value_at(2*kx  ,ky,  kx+ky*nx/2);
    tiBFFT.set_value_at(2*kx+1,ky,-(kx+ky*nx/2));
# Now the real and complex values are negatives
#     of one another
# the minimum magnitude is zero
# the maximum magnitude is 

#  This does not correspond to a
#  real image, but that is no matter for this section



# Let's retrieve the information using get_complex_at
# kx =0, ky=0, kx+4*ky=0
tiBFFT.get_complex_at(0,0)  # gives 0 correctly

# kx =1, ky=1,, kx+4*ky=5
tiBFFT.get_complex_at(1,1);   #  gives 4 -4j correctly
tiBFFT.get_value_at(2,1); tiBFFT.get_value_at(3,1);  # is 4, -4

# kx =3, ky=-1 ( or 5),, kx+4 ky=23
tiBFFT.get_complex_at(2,4-5)  # should be 14 - 14j
tiBFFT.get_value_at(4,4); tiBFFT.get_value_at(5,4);  # is 14, -14


#  Friedel pairs
# kx =-1, ky=-1,
#  This is Friedel pair of kx=1,ky=1
# kx=1 , ky=1,   kx+3*ky=4
tiBFFT.get_complex_at(-1,-1);   #  gives 4 + 4j correctly  
tiBFFT.get_value_at(2,1); tiBFFT.get_value_at(3,1);  # is 4, -4


###########################################
##  Part 2.   Counting numbers of variables
#     get_complex_at   for N odd, 3D

# Naively there are a total of 150= (N+1)(N)N real values in the complex array
#  that stores the 6 by 5 by 5 FFT image
# However the original image housed only   5*5*5 =125 real values


# The 25= N^2 redundancies are 
# kx=0 , ky= 0 , kz= 0   value is  real        1

# [kx=0, ky=0]     kz=1  or kz=4  and kz=2 kz=3    Friedl pairs         2*2
# [kx=0, kz=0]     ky=1  or ky=4  and ky=2 ky=3    Friedl pairs         2*2
# [ kx =0, and (ky,kz) = (1:4,1:4) ]  has  (N-1)^2=16 entries
#                      half are Friedel pairs of the other half =>     16

# Generally there are  1 reals
#  Friedel pair redundancies    2  N of the type of the first two lines of last paragraph
#  Friedel pairs redundancies   (N-1)*(N-1)   of the last four lines of last paragraph

# This gives N^2 = 25 conditions




N=5
ti=EMData(N,N); ti.to_zero();
for jx in range(N):
  for jy in range(N):
    ti.set_value_at(jx,jy, random());

tiB    = ti.process('xform.phaseorigin.tocorner')
tiBFFT = tiB.do_fft();


#  1 special values
print tiBFFT.get_complex_at(0,0);     #  yields real value


#  (N -1)  =4 redundancies
kx=0;
for ky in range(1,(N+1)/2): #  (N/2-1) values
    print tiBFFT.get_complex_at(kx, ky), tiBFFT.get_complex_at(kx,-ky);
    print 'Last two values should be complex conjugates';

# so real values = complex values
# N^2= (N+1)(N) -1   -((N-1)/2)(2) 

###########################################################
#   Part 3  set_complex_at
#
#  So there are  N special points to worry about when using set_complex_at
#  1  of them we must make sure are real
#  N-1 of them we must set the Friedel pair

#   This is tricky, because  the set command may influence
#    multiple Points
#  So if we set anything with kx=0   
#      we have to ensure that something else is true
#   in order that our complex image is an FFT of a real image
#  (the Friedel pair must be a complex conjugate)

# set_complex_at will accept those parts of the set that are admissible
#  For examples

tiBFFT.set_complex_at(2,0,1+2j)
tiBFFT.get_complex_at(2,0)#  gives 1

# set_complex_at won't do wrapping 
tiBFFT.set_complex_at(2,3,1+2j)

tiBFFT.get_complex_at(2,3)# Out[62]: 0j
tiBFFT.get_complex_at(2,-2) # will be unaffected by the last set

# the redundant information it will create correctly

tiBFFT.set_complex_at(0,1,1+2j)
tiBFFT.get_complex_at(0,1)#   1 +2j
tiBFFT.get_complex_at(0,-1)#  1-2j


# Let's check this completely

N=5;
kx=0  
for ky in range(1,(N+1)/2):#  (N/2-1) values
    Value =complex(1,ky);
    tiBFFT.set_complex_at(kx,ky,Value);
    print tiBFFT.get_complex_at(kx,ky),tiBFFT.get_complex_at(kx,-ky);# this is the cc of above line
    print 'Last two values should be complex conjugates';




