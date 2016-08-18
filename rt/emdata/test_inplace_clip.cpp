
#include "emdata.h"
#include "geometry.h"
using namespace EMAN;

#include <vector>
using std::vector;

#include <iostream>
using std::cout;
using std::endl;


void test_read_images()
{
	EMData* data = new EMData;
	vector<int> idxs;
	vector<EMData* > v = data->read_images("test_out1.img", idxs, false);

	
	cout << "v had " << v.size() << " images" << endl;
	
	for ( vector<EMData* >::const_iterator it = v.begin(); it != v.end(); ++it)
	{
		(*it)->write_image("test_out.img", -1 );
	}
}

#include<ctime>


void test_clip_inplace_time()
{
	int n = 6;

	Region larger_region(-128, -128, -128, 768, 768, 768);
	Region smaller_region(128, 128, 128, 256, 256, 256);
	clock_t t1=clock();
	for( int i = 0; i < n; ++i )
	{		
		EMData* data = new EMData("threed_large.mrc");
		data->clip_inplace(larger_region);
		delete data;
	}
	clock_t t2=clock();

	cout<<"Inplace larger clips took " << double(t2-t1)/CLOCKS_PER_SEC << " seconds" << endl;

	t1=clock();
	for( int i = 0; i < n; ++i )
	{		
		EMData* data = new EMData("threed_large.mrc");
		EMData* clipped = data->get_clip(larger_region);
		delete data;
		delete clipped;
	}
	t2=clock();

	cout<<"Out of place larger clips took " << double(t2-t1)/CLOCKS_PER_SEC << " seconds" << endl;

	t1=clock();
	for( int i = 0; i < n; ++i )
	{		
		EMData* data = new EMData("threed_large.mrc");
		data->clip_inplace(smaller_region);
		delete data;
	}
	t2=clock();

	cout<<"Inplace smaller clips took " << double(t2-t1)/CLOCKS_PER_SEC << " seconds" << endl;

	t1=clock();
	for( int i = 0; i < n; ++i )
	{		
		EMData* data = new EMData("threed_large.mrc");
		EMData* clipped = data->get_clip(smaller_region);
		delete data;
		delete clipped;
	}
	t2=clock();

	cout<<"Out of place smaller clips took " << double(t2-t1)/CLOCKS_PER_SEC << " seconds" << endl;
}

void test_clip_inplace_visually()
{

	Region region(-10, 30, 23, 116, 133, 110);
	EMData* data = new EMData("threed.mrc");
	data->clip_inplace(region);
	data->write_image("clipped_inplace.mrc");
}

bool test_clip_inplace( const Region& region)
{
	EMData* data = new EMData("threed_small.mrc");

	EMData* clipped = data->get_clip(region);
	data->clip_inplace(region);
	EMData* tmp = ((*clipped) - (*data));
	int z = tmp->get_zsize();
	int y = tmp->get_ysize();
	int x = tmp->get_xsize();
	//cout << "The output dimensions are " << x << " " << y << " " << z;
	bool failure = false;
	for( int i = 0; i < z; ++i )
		for ( int j = 0; j < y; ++j)
			for ( int k = 0; k < x; ++k )
			{
				float score = tmp->get_value_at(k,j,i);
				if ( score != 0 )
				{
					failure = true;
					cout << "Error, the pixel value was not zero, it was " << score << " at [x,y,z] " << x << " " << y << " " << z << endl;
				}
			}

	

	delete data;
	delete clipped;
	delete tmp;

	return failure;
}

void test_many_clip_inplace()
{
	vector<Region> bad_regions;
	for ( int x0 = -1; x0 <= 1; x0 +=1 )
	for ( int y0 = -1; y0 <= 1; y0 +=1 )
	for ( int z0 = -1; z0 <= 1; z0 +=1 )
	for ( int xsize = 31; xsize <= 33; xsize +=1 )
	for ( int ysize = 31; ysize <= 33; ysize +=1 )
	for ( int zsize = 31; zsize <= 33; zsize +=1 )
	{
		cout << "Inplace testing with translataion [" << x0 << "," << y0 << "," << z0 << "] and dimensions [" << xsize << "," << ysize << "," << zsize << "]";
		Region region(x0,y0,z0,xsize,ysize,zsize);
		if ( test_clip_inplace(region) )
		{
			cout << "....................FAIL" << endl;
			bad_regions.push_back(region);
		}
		else 
		{
			cout << "....................ok" << endl;
		}
	}
	
	cout << endl << "Overall there were " << bad_regions.size() << " failures" << endl;
}

int main(int argc, char** argv)
{

	test_many_clip_inplace();
	//test_clip_inplace_visually();
// 	test_clip_inplace_time();

	return 1;
}
