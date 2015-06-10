/*
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by.edge
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
 */

#include "tomoseg.h"
#include <algorithm>
using namespace EMAN;

float DistToLine(int x0,int y0,int x1,int y1,int x2,int y2){
	return abs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))/sqrt((float)((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
}

TomoObject::TomoObject(vector<Vec3i> allpt, float md, int slice){
	nowslice=slice;
	maxdist=md;
	// find one end point
	for (int i=0; i<(int)allpt.size(); i++){
		if (allpt[i][2]==1){
			allpoints.push_back(allpt[i]);
			allpt[i][2]=-1;
			break;
		}
	}
	
	// sort points from one endpoint to the other
	for (int i=0; i<(int)allpt.size()-1; i++){
		float mindst=1e10;
		int nxti=-1;
		Vec3i nowp=allpoints.back();
		for (int j=0; j<(int)allpt.size(); j++){
			if (allpt[j][2]>0){
				float dst=(allpt[j][0]-nowp[0])*(allpt[j][0]-nowp[0])+(allpt[j][1]-nowp[1])*(allpt[j][1]-nowp[1]);
				if (dst<mindst){
					mindst=dst;
					nxti=j;
				}
			}
		}
		if (nxti>=0){
			allpoints.push_back(allpt[nxti]);
			allpt[nxti][2]=-1;
		}
		else
			break;
	}
// 	for (int i=0; i<(int)allpoints.size(); i++)
// 		printf("(%d,%d)",allpoints[i][0],allpoints[i][1]);
// 	printf("\n");
	ptid.push_back(0);
	ptid.push_back((int)allpoints.size()-1);
	
	
	// break into line segments
	while(1){
		int newpt=bend();
		if (newpt>0){
			ptid.push_back(newpt);
			sort(ptid.begin(),ptid.end());
// 			for (int i=0; i<ptid.size(); i++) printf("%d\t",ptid[i]);
// 			printf("\n");
			
		}
		else
			break;
	}
	
}

int TomoObject::bend(){
	// distance to current line segments
	float maxd=maxdist;
	int newpt=-1;
	for (int i=1; i<(int)ptid.size(); i++){
		int x1=allpoints[ptid[i-1]][0],y1=allpoints[ptid[i-1]][1];
		int x2=allpoints[ptid[i]][0],y2=allpoints[ptid[i]][1];
		for (int j=ptid[i-1]; j<ptid[i]; j++){
			int x0=allpoints[j][0],y0=allpoints[j][1];
			float dst=DistToLine(x0,y0,x1,y1,x2,y2);
			if (dst>maxd){
				maxd=dst;
				newpt=j;
			}
		}
		
	}
	return newpt;
}

void TomoObject::write_imod(FILE *fp){
	fprintf(fp,"color 0 1 0 0\nopen\ncontour 0 0 %d\n",(int)ptid.size());
	for (int i=0; i<(int)ptid.size(); i++){
		fprintf(fp,"%d %d %d\n",allpoints[ptid[i]][0],allpoints[ptid[i]][1],nowslice);
	}
}

int TomoObject::get_size(){
	return (int)ptid.size();
}

void TomoObject::enlong(EMData *bwmap,EMData *skelmap){
	// first push the current segment in.
	int myid=skelmap->get_value_at(allpoints[0][0],allpoints[0][1]);
	segid.push_back(myid);
	
	//
	
}

int TomoSeg::read_skelmap(EMData *map){
	if (!map)
		return 0;
	skelmap=map;
	return 1;
}

int TomoSeg::generate_objects(int numo, float maxdist, int nowslice){
	if (!skelmap){
		printf("No map input.\n");
		return 0;
	}

	// Calculate area of each segment	
	EMData *rank=skelmap->process("threshold.notzero");
	rank->process_inplace("morph.object.density",Dict("thresh",0,"more_neighbor",true));
	int nx=rank->get_xsize();
	int ny=rank->get_ysize();
	
	// Take numo largest objects
	for (int i=0; i<numo; i++){
		float max=rank->get_attr("maximum");
		vector <Vec3i> pts;
		float curid=-1;
		bool branch=false;
		if (max<=0)
			break;
		for (int x=0; x<nx; x++){
			for (int y=0; y<ny; y++){
				if (rank->get_value_at(x,y)==max){
					if (curid<0)
						curid=skelmap->get_value_at(x,y);
					if (curid==skelmap->get_value_at(x,y)){
						int nb=check_neighbors(x,y);
						if (nb>2) branch=true;
						pts.push_back(Vec3i(x,y,nb));
						rank->set_value_at_fast(x,y,0);
					}
				}
			}
		}
		if(verb) printf("id=%d,  size=%d",int(curid),int(pts.size()));
		if (branch){
			if (verb) printf("\n\thave branch, throw out.\n");
		}
		else{
			objs.push_back(TomoObject(pts,maxdist,nowslice));
			if (verb) printf("   \t%d points\n",objs.back().get_size() );
		}
		
	}
	
	delete rank;
	return 1;
}

void TomoSeg::write_imod(const char *file){
	FILE *fp = fopen(file, "w");
	int nx=skelmap->get_xsize();
	int ny=skelmap->get_ysize();
	int nz=skelmap->get_zsize();
	fprintf(fp,"imod %d\n",(int)objs.size());
	fprintf(fp,"max %d %d %d\n",nx,ny,nz);
	for (int i=0; i<(int)objs.size(); i++){
		fprintf(fp,"object %d 1 0\n",i);
		objs[i].write_imod(fp);
		fprintf(fp,"\n");
	}
	fclose(fp);
}

int TomoSeg::check_neighbors(int x,int y){
	int nb=-1;
	for (int i=-1; i<=1; i++){
		for (int j=-1; j<=1; j++){
			nb+= (skelmap->get_value_at(x+i,y+j)>0 ? 1 : 0);
		}
	}
	return nb;
}
