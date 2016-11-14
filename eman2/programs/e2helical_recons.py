#!/usr/bin/env python
from EMAN2 import *
import random
from math import *
import os
from os import system
import sys
from e2simmx import cmponetomany
from reconstruction import recons3d_4nn
import numpy # import numpy namespace

def main():
  progname=os.path.basename(sys.argv[0])
  usage="""prog [options] blablabla"""
  parser = EMArgumentParser(usage=usage,version=EMANVERSION)
  parser.add_argument("--input",dest="input",default=None,type=str,help="particles data") #need this input as particles set 
  parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
  parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness") 
  parser.add_argument("--parallel",type=str,help="Parallelism string",default=None)
 
  global options
  (options,args)=parser.parse_args()
 
############below are parameter sets########## 

  proj0='bdb:projections#proj_02142011' # this is a string object, is the output projections set database name
  simx1='bdb:simmx#simmx_02142011' # this is simmilarity matrix file name
  #aligned1='bdb:aligned_particle#aligned_ptcl_02112011'
  update_vol1='bdb:updated_models#new_model_02142011'
  least_output1="least_data02142011.dat"
  
############above are paramter sets#######################  
  ptcls=EMData.read_images(options.input) # read in the original particle set, this typically does not change 
  print options.input  # input is the original particles set
  print "total %d particles"%len(ptcls) # print out how many particles in the set
  reso=90 # define projections number, more projections mean higher resolution
  #ret=EMData.read_images(update_vol1) # read in the volume set
  orts=[Transform({"type":"eman","az":360/reso*j,"phi":0,"alt":90}) for j in xrange(reso)]
   
  for cycle1 in xrange(100):   # the model refinement cycle
  #### the projections has to be cleared before next run
    
    cycle=len(ret)+2*cycle1
    print "this is cycle %d" %cycle
    proj1=proj0+"cycle_%d" %cycle
    print proj1
    ret=EMData.read_images(update_vol1) # volume set is updated, need reload
    print "len(ret)=%d" %len(ret)
    
    projs=[ret[cycle].project("standard",ort) for ort in orts] # orts is the list of transformations, projs is the projections list corresponding those transformations
    
    try: os.mkdir("projections")
    except: pass
    #dcts=db_list_dicts("bdb:projections")
    for i in xrange(reso):
      print "this is projection %d" %i

      projs[i].write_image("%s" %proj1,-1) # write out the 90 projections of volume to file for the following e2simmx.py to use
    #os.system("rm -rf simmx") # clean the simmx folder
    
    try: os.mkdir("simmx")
    except: pass
    
    launch_childprocess("e2simmx.py %s  %s %s -f --saveali --cmp=ccc --align=rotate_translate_flip --aligncmp=ccc --verbose=2  --parallel=thread:4" %(proj1, options.input, simx1)) 
  # e2simmx.py compare the 90 projections with the particles set to get a similarity matrix to store the similarity score and corresponding transformation parametes xform.projection for the back projection program recons3d_4nn() to use
  ###################################  
    
    score=EMData.read_images("%s" %simx1) # read the scrore from similarity matrix simmx_02022011a 
    bslst=[]
    for i in xrange(len(ptcls)):
      bs=[]
      bs=[score[0][j] for j in xrange(i*reso,(i+1)*reso)] # for each particle find all scores corresponding to each projection, the score is stored in the first line score[0][j] (score[0][j] means the score of praticle i and jth projection).
      best_score=min(bs) # bs is the list of scores of ith particle and the 90 projections, best_score is the smallest number of the scrore in bs  
      bslst.append((best_score,i)) # bslst is the best score list of all the particles, ith particle has a best score and stored in bslst and so on for the other particle.
      n=bs.index(best_score) # find the projection number n corresponding best score of particle i.
      ptcls[i]["match_n"]=n # assign n to particle i
      if i%30==0: print "%d, %2d%% complete!" %(i,i*100/len(ptcls))
      ptcls[i]["match_qual"]=best_score 
      ptcls[i]["xform.projection"]=orts[n]

    #try: os.mkdir("aligned_particle") #created a folder to put the aligned particles
    #except: pass
    bslst.sort() # sort the best score list from smallest to bigger (-1 is best)
    #bslst.reverse()
    aptcls=[] # aligned particles set
    for i in xrange(int(len(ptcls)*60/100.0)): # use the 70% total particles which has better best_score
      n=ptcls[bslst[i][1]]["match_n"] # the ith best score corresponds particle bslst[i][1]
      aptcls.append(ptcls[bslst[i][1]].align("rotate_translate_flip",projs[n],{},"dot",{}))
      #aptcls[-1].process_inplace("xform.centerofmass",{})
      if i%30==0: print "score=%f, %d, %2d%% complete!" %(bslst[i][0], i,i*100*8/6/len(ptcls))
      #aptcls[i].process_inplace("xform.centerofmass",{})
      #aptcls[i].write_image(aligned1,-1) # write out the aligned particle set

    vol=recons3d_4nn(aptcls) # recons3d_4nn(aptcls) use the aligned particles set aptcls to generate a volume using back projection
    try: os.mkdir("updated_models") 
    except: pass
    vol.write_image(update_vol1,-1) # write out the volume file to the end of list
    print "back projection completed"
   # display(vol) 
    
    new_volume=vol.copy()  # new_volume is the volume in the volume set new_vol
    nx=new_volume.get_xsize() # the boxsize of the new volume new_vol[0]  
############################################
######################################################### 
    dz=15.06 # initial guessing of dz
    phi_array=[]
    std_array=[]
    for ite in xrange(1): 
      print "this is %d cycle" %ite
      d_phi=35.+0.05*ite # initial guessing of d_phi, the unit is degree
  ################## change the volume to cyliner
      '''
      for i in xrange(nx):
        for j in xrange(nx):
          for k in xrange(nx):
            if sqrt((i-nx/2-0.5)*(i-nx/2-0.5)+(j-nx/2-0.5)*(j-nx/2-0.5))>nx/2:
              new_volume.set_value_at(i,j,k,0.0) # assign value 0.0 to the voxel (i,j,k) for the distance>nx/2
      '''
      avg_vol=EMData(nx,nx,nx) # initialize the averaged volume for 9 subunits
      mean_std=0.0 # intialize mean squared deviation
      for k in xrange(int (nx-dz)/2, int (nx+dz)/2): # k ranged from lower plane to upper plane ---> width of central volume
        for i in xrange(nx):
          for j in xrange(nx):
            if new_volume.get(i,j,k)!=0.0: # the voxel should not be 0.0
              sublst=[] # sublst is the list containning the 9 voxel value of the neighboring subunits
              r=sqrt((i-nx/2-0.5)*(i-nx/2-0.5)+(j-nx/2-0.5)*(j-nx/2-0.5))  # find r, the distance to the center (nx/2-0.5,nx/2-0.5)
              theta=180.0*atan((j-nx/2-0.5)/(i-nx/2-0.5))/pi # find the angle theta ranging from 0 to 360
              if i<nx/2:
	        theta+=180.0 # this is very important, since atan(theta) range from -pi/2 to pi/2
              #print "r=%f\t theta=%f" %(r,theta) 
              for units in xrange(-4,5): #use 8 neighboring subunits from -4 to +4: -4 -3 -2 -1 0 1 2 3 4
                if units==0:
                  sublst.append(new_volume.get(i,j,k)) # put the central subunit to the list      
                else:
                  temp_theta=(theta+units*d_phi)*pi/180.0
                  temp_x=int(r*cos(temp_theta)+nx/2+0.5)
                  temp_y=int(r*sin(temp_theta)+nx/2+0.5)
                  temp_z=int(k+units*dz)
                  sublst.append(new_volume.get(temp_x,temp_y,temp_z))
           
              temp_mean=numpy.mean(sublst) # mean value
              avg_vol.set_value_at(i,j,k,temp_mean)
              temp_std=numpy.std(sublst)
              mean_std+=temp_std 
        print "i,j,k=%d,%d,%d\t std=%f" %(i,j,k,mean_std)
      print "dz=%f\t d_phi=%f\t mean_deviation=%f\n" %(dz,d_phi,mean_std)
    
      thisfile=file(least_output1,'a') # output the least square data to file
      print>>thisfile,dz,d_phi,mean_std
      thisfile.close()
      phi_array.append(d_phi)
      std_array.append(mean_std)

    results=numpy.polyfit(phi_array,std_array,2) # use numpy.polyfit to get least square fit of phi and std
    fine_phi=-0.5*results[1]/results[2]
   # print "fine_phi=%f" %fine_phi
###############################################
  
    xf=Transform() # needed to run get_sym
    xf.to_identity()
   # display(avg_vol)
    avg_vol.translate(0,0,-(nx-dz)/2+1) # translate avg_vol to the bottom just 1 voxel above
    #display(avg_vol)
    dcopy=avg_vol.copy() # dcopy is the copy of the averaged volume
    sym_str="h4,20,%d,%d" %(d_phi,dz) # define h-symmetry parameter string
    for i in xrange(17):
      temp_av=avg_vol.copy()
      temp_av.transform(xf.get_sym(sym_str,i)) # dc is the h-symmetrized subunit
      dcopy.add(temp_av)
    dcopy.sub(avg_vol)  
    dcopy.mult(1.0/4.0) # cyclic 4 symmetry, density multiplied 4, should change back
    dcopy.translate(0,0,-2) #shift downward 2 pixels
##################### mask the model ###################
    '''
    for i in xrange(nx):
      for j in xrange(nx):
        for k in xrange(nx):
          #if sqrt((i-nx/2-0.5)*(i-nx/2-0.5)+(j-nx/2-0.5)*(j-nx/2-0.5))>nx/2*0.9:
          if dcopy.get(i,j,k)==0.0: 
            dcopy.set_value_at(i,j,k,0.000001) # assign value 0.0 to the voxel (i,j,k) for the distance>nx/2
    '''
############ end of mask ##############  
    #display(dcopy) 
    try: os.mkdir("updated_models")
    except: pass
    dcopy.write_image(update_vol1,-1)

if __name__=="__main__":  
  main()

