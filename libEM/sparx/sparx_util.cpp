#include "sparx_util.h"
#include <list>
#include <vector>

using namespace EMAN;
using std::list;
using std::vector;

bool SparxUtil::peakcmp(const Pixel& p1, const Pixel& p2) {
    return (p1.value > p2.value);
}

ostream& operator<< (ostream& os, const Pixel& peak) {
    os <<  peak.x <<  peak.y << peak.z
       << peak.value;
    return os;
}

vector<float> SparxUtil::peak_search(EMData* img, int ml, float invert)
{
 	 EMData& buf = *img;
	 vector<Pixel> peaks;
 	 int img_dim;
 	 int i,j,k,itx,ity,itz;
 	 int i__1,i__2;
 	 int j__1,j__2;
 	 int k__1,k__2;
 	 bool peak_check;
 	 img_dim=buf.get_ndim();
 	 vector<float>ix,jy,kz,res;
 	 int nx = buf.get_xsize();
 	 int ny = buf.get_ysize();
 	 int nz = buf.get_zsize();
 	 switch (img_dim)
     {
 	 case(1):
		for(i=0;i<=nx-1;++i)
 	  	{
 	   		i__1=(i-1+nx)%nx;
 	   		i__2=(i+1)%nx;
 	      	peak_check=buf(i)*invert>buf(i__1)*invert && buf(i)*invert>buf(i__2)*invert;
	 	  	if(peak_check)
		  		{peaks.push_back(Pixel(i, 0, 0, buf(i)*invert));}  
	 	}
 	 break;
 	 case(2):
		for(j=0;j<=ny-1;++j)
 	    {  
 	    	j__1=(j-1+ny)%ny;
 		 	j__2=(j+1)%ny;
 	        for(i=0;i<=nx-1;++i)
			{ 
				i__1=(i-1+nx)%nx;
			  	i__2=(i+1)%nx;
			  	peak_check=(buf(i,j)*invert>buf(i,j__1)*invert) && (buf(i,j)*invert>buf(i,j__2)*invert);
			  	peak_check=peak_check && (buf(i,j)*invert>buf(i__1,j)*invert) && (buf(i,j)*invert>buf(i__2,j)*invert);
			  	peak_check=peak_check && (buf(i,j)*invert>buf(i__1,j__1)*invert) && ((buf(i,j)*invert)> buf(i__1,j__2)*invert);
			  	peak_check=peak_check && (buf(i,j)*invert>buf(i__2,j__1)*invert) && (buf(i,j)*invert> buf(i__2,j__2)*invert);
			  	if(peak_check)
 			    {peaks.push_back(Pixel(i, j, 0, buf(i,j)*invert));}
			}
 		}
 	 break;
 	 case(3):
		for(k=0;k<=nz-1;++k)
	   	{  
	   		kz.clear();
	      	k__1=(k-1+nz)%nz;
	      	k__2=(k+1)%nz;
	      	kz.push_back(k__1);
	      	kz.push_back(k);
	      	kz.push_back(k__2);
 	     	for(j=0;j<=ny-1;++j)
		    {
		    	jy.clear();
		     	j__1=(j-1+ny)%ny;
 	            j__2=(j+1)%ny;
		     	jy.push_back(j__1);
		     	jy.push_back(j);
		     	jy.push_back(j__2);
 	            for(i=0;i<=nx-1;++i)
			   	{ 
			   		ix.clear();
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
						peaks.push_back(Pixel(i, j, k, buf(i,j,k)*invert));
					} 
			      }
			   }
			}
		break;
     	}
   	     	    	 
   sort(peaks.begin(),peaks.end(), peakcmp);   
   int count=0;
   float xval=(*peaks.begin()).value;
  for (vector<Pixel>::iterator it = peaks.begin(); it != peaks.end(); it++)  
  {count=count+1;
       if(count<=ml)
	  {
	  res.push_back((*it).value);
	  res.push_back((*it).x);
	  if(img_dim!=1)
	       {res.push_back((*it).y);}
	  if(nz!=1)
		{res.push_back((*it).z);}
	  res.push_back((*it).value/xval);
	  res.push_back((*it).x-float(nx)/2.f);
	  if(img_dim!=1)
	      {res.push_back((*it).y-float(ny)/2.f);}
	  if(nz!=1)
	     {res.push_back((*it).z-float(nz)/2.f);} 
	  
	}
   }  
  res.insert(res.begin(),1,img_dim);
  return res;
 }	
	  

