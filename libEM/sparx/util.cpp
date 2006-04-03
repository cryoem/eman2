#include "emdata.h"
using namespace EMAN;
namespace EMAN {
    
    EMData*    
    struct Peak {
    float val;
    float xpos;
    float ypos;
    float zpos;
    Peak (float val_, float xpos_=0, float ypos_=0, float zpos_=0)
        : val(val_), xpos(xpos_), ypos(ypos_), zpos(zpos_) {}
};

bool peakcmp(const Peak p1, const Peak p2) {
    return (p1.val > p2.val);
}
ostream& operator<< (ostream& os, Peak peak) {
    os <<  peak.xpos <<  peak.ypos << peak.zpos
       << peak.val;
    return os;
}
     vector<float> EMData::peak_search(int ml, float invert)
{
 	 EMData& buf = *this;
	 list<Peak> peaks;
 	 int img_dim;
 	 int i,j,k,itx,ity,itz;
 	 int i__1,i__2;
 	 int j__1,j__2;
 	 int k__1,k__2;
 	 bool peak_check;
 	 img_dim=get_ndim();
 	 vector<float>ix,jy,kz,res;
 	 switch (img_dim)
      {
 	 case(1):
 	   for(i=0;i<=nx-1;++i)
 	  {
 	   i__1=(i-1+nx)%nx;
 	   i__2=(i+1)%nx;
 	      peak_check=buf(i)*invert>buf(i__1)*invert && buf(i)*invert>buf(i__2)*invert;
	 	  if(peak_check)
		  {peaks.push_back(Peak(buf(i)*invert,i));}  
	 }
 	 break;
 	 case(2):
 	   for(j=0;j<=ny-1;++j)
 	      {  j__1=(j-1+ny)%ny;
 		 j__2=(j+1)%ny;
 	           for(i=0;i<=nx-1;++i)
			{ i__1=(i-1+nx)%nx;
			  i__2=(i+1)%nx;
			  peak_check=(buf(i,j)*invert>buf(i,j__1)*invert) && (buf(i,j)*invert>buf(i,j__2)*invert);
			  peak_check=peak_check && (buf(i,j)*invert>buf(i__1,j)*invert) && (buf(i,j)*invert>buf(i__2,j)*invert);
			  peak_check=peak_check && (buf(i,j)*invert>buf(i__1,j__1)*invert) && ((buf(i,j)*invert)> buf(i__1,j__2)*invert);
			  peak_check=peak_check && (buf(i,j)*invert>buf(i__2,j__1)*invert) && (buf(i,j)*invert> buf(i__2,j__2)*invert);
			  if(peak_check)
 			      {peaks.push_back(Peak(buf(i,j)*invert,i,j));}
			}
 		}
 	 break;
 	 case(3):
         for(k=0;k<=nz-1;++k)
	   {  kz.clear();
	      k__1=(k-1+nz)%nz;
	      k__2=(k+1)%nz;
	      kz.push_back(k__1);
	      kz.push_back(k);
	      kz.push_back(k__2);
 	     for(j=0;j<=ny-1;++j)
		    {jy.clear();
		     j__1=(j-1+ny)%ny;
 	             j__2=(j+1)%ny;
		     jy.push_back(j__1);
		     jy.push_back(j);
		     jy.push_back(j__2);
 	              for(i=0;i<=nx-1;++i)
			   { ix.clear();
			     i__1=(i-1+nx)%nx;
			     i__2=(i+1)%nx;
			     ix.push_back(i__1);
			     ix.push_back(i);
			     ix.push_back(i__2);
			     peak_check=true ;
				 for(itx=0;itx<= 2; ++itx)
				     {  
				     
				         for(ity=0;ity<= 2; ++ity)
					    {  
					    
					         for(itz=0;itz<= 2; ++itz) 
						 
						 { 
			if((buf(i,j,k)*invert<=buf(ix[itx],jy[ity],kz[itz])*invert) && i !=ix[itx] && j !=jy[ity] && k !=kz[itz])			   	
						        {peak_check=false;}						   
						 
						  }
	
					     }		 
			            }
					    if(peak_check)
				          {
					   peaks.push_back(Peak(buf(i,j,k)*invert,float(i),float(j),float(k)));
					  } 
					
			      }
				  
			   }
		}	     	    	 
 	 break;
	 }
}
