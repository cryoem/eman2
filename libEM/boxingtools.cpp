/**
 * $Id$
 */

/*
 * Author: David Woolford, 04/15/2008 (woolford@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 *
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
 * it under the terms of the GNU General Public License as published by
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
 *
 * */


#include "boxingtools.h"
#include "exception.h"
using namespace EMAN;


#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

// For random
#include <cstdlib>
#include <ctime>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <string>
using std::string;

#include <algorithm>
// find, min_element

vector<Vec3f> BoxingTools::colors = vector<Vec3f>(); // static init
BoxingTools::CmpMode BoxingTools::mode = SWARM_RATIO;

#define SVD_CLASSIFIER_DEBUG 0


#if SVD_CLASSIFIER_DEBUG
 void printMatrix( const gsl_matrix * const A, const unsigned int rows, const unsigned int cols, const string& title = "" )
{
	cout << "Printing matrix " << title << endl;
	cout << "It has " << rows << " rows and " << cols << " columns " << endl;
	for( unsigned int i = 0; i < rows; ++i )
	{
		for( unsigned int j = 0; j < cols; ++j )
		{
			cout << gsl_matrix_get(A,i,j) << " ";
		}
		cout << endl;
	}

}

void printVector( const gsl_vector * const A, const unsigned int length, const string& title = "" )
{
	cout << "Printing vector " << title << endl;
	for( unsigned int i = 0; i < length; ++i )
	{
		cout << gsl_vector_get(A,i) << " ";
	}
	cout << endl;
}

void print_map( const map<unsigned int, unsigned int>& mapping )
{
	for(map<unsigned int, unsigned int>::const_iterator it = mapping.begin(); it != mapping.end(); ++it)
	{
		cout << it->first << " " << it->second << endl;
	}
}
#endif

BoxSVDClassifier::BoxSVDClassifier(const vector<vector<float> >& data, const unsigned int& classes) :
		mData(data), mClasses(classes)
{
	setDims( mData );
}


BoxSVDClassifier::~BoxSVDClassifier()
{

}

bool BoxSVDClassifier::setDims( const vector<vector<float> >& data )
{
	mColumns = mData.size();
	vector<vector<float> >::const_iterator it = data.begin();
	mRows = it->size();
	it++;
	for( ; it != data.end(); ++it )
	{
		if ( it->size() != mRows )
		{
			cerr << "ERROR: can not initial the BoxSVDClassifier with vectors of un-equal lengths " << endl;
			cerr << "The vector lengths that did not agree were " <<  mRows << " and " << it->size() << endl;
			return false;
		}
	}

	return true;
}


map< unsigned int, unsigned int> BoxSVDClassifier::go()
{
	//	This is done in the constructor
	// 	setDims(mData);


	unsigned int local_columns = mColumns;
	if ( mRows < mColumns )
	{
// 		cerr << "Warning: gsl SVD works only when m > n, you have m = " << mRows << " and n = " << mColumns << endl;
		// This local adaptation means things will proceed the same way even if there are more columns in A then rows
		// Every input data is still classified, just the SVD eigenvectors are found using a subset of all the data
		local_columns = mRows;
	}

	gsl_matrix * U = gsl_matrix_calloc( mRows, local_columns );
	gsl_matrix * A = gsl_matrix_calloc( mRows, mColumns );
	for ( unsigned int i = 0; i < mRows; ++i )
	{
		for ( unsigned int j = 0; j < mColumns; ++j )
		{
			gsl_matrix_set( A, i, j, mData[j][i] );
			if ( j < local_columns )
				gsl_matrix_set( U, i, j, mData[j][i] );
		}
	}
#if SVD_CLASSIFIER_DEBUG
	printMatrix( A, mRows, mColumns, "A" );
#endif

	gsl_matrix * V = gsl_matrix_calloc( local_columns, local_columns );
	gsl_vector * S = gsl_vector_calloc( local_columns );
	gsl_vector * work = gsl_vector_calloc( local_columns );

	if ( gsl_linalg_SV_decomp (U, V, S, work) )
	{
		cerr << "ERROR: gsl returned a non zero value on application of the SVD" << endl;
	}

#if SVD_CLASSIFIER_DEBUG
	printMatrix( U, mRows, local_columns, "U" );
	printVector( S, local_columns, "S" );
	printMatrix( V, local_columns, local_columns, "V");
#endif

	// normalize the columns of matrix A
	for ( unsigned int j = 0; j < mColumns; ++j )
	{
		float norm = 0;
		for ( unsigned int i = 0; i < mRows; ++i )
		{
			norm += (float)(gsl_matrix_get( A, i, j)*gsl_matrix_get( A, i, j));
		}
		norm = sqrtf(norm);
		for ( unsigned int i = 0; i < mRows; ++i )
		{
			gsl_matrix_set( A, i, j, gsl_matrix_get(A,i,j)/norm);
		}
	}

#if SVD_CLASSIFIER_DEBUG
	for ( unsigned int j = 0; j < mColumns; ++j )
	{
		double norm = 0;
		for ( unsigned int i = 0; i < mRows; ++i )
		{
			norm += gsl_matrix_get( A, i, j)*gsl_matrix_get( A, i, j);
		}
		cout << "For column " << j << " the squared norm is " << norm << endl;
	}
#endif


	gsl_matrix * svd_coords = gsl_matrix_calloc( mColumns, mColumns );
	// Correlate the columns of A with the columns of U and store the information in a martrix called svd_coords
	for ( unsigned int i = 0; i < mColumns; ++i )
	{
		for ( unsigned int j = 0; j < local_columns; ++j )
		{
			double result = 0.0;
			for ( unsigned int k = 0; k < mRows; ++k )
			{
				result += gsl_matrix_get(A,k,i)*gsl_matrix_get(U,k,j);
			}
			gsl_matrix_set( svd_coords, i, j, result);
		}
	}

#if SVD_CLASSIFIER_DEBUG
	printMatrix( svd_coords, mColumns, mColumns, "svd_coords" );
#endif

	map< unsigned int, unsigned int> grouping = randomSeedCluster(svd_coords, mColumns);

	for ( unsigned int i = 0; i < 20; ++ i )
	{
		grouping = getIterativeCluster(svd_coords, grouping);
	}

	gsl_matrix_free(A);
	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_matrix_free(svd_coords);

	return grouping;
}

map< unsigned int, unsigned int> BoxSVDClassifier::getIterativeCluster(const gsl_matrix* const svd_coords, const map< unsigned int, unsigned int>& current_grouping)
{
	// Space to store the reference vectors
	gsl_matrix * ref_coords = gsl_matrix_calloc( mClasses, mColumns );

	// Assumes there are a total of mClasses in the current_groupings mapping
	for(unsigned int i = 0; i < mClasses; ++i)
	{
		unsigned int tally = 0;
		for (map< unsigned int, unsigned int>::const_iterator it = current_grouping.begin(); it != current_grouping.end(); ++it )
		{
			if ( it->second == i )
			{
				for( unsigned int j = 0; j < mColumns; ++j )
				{
					gsl_matrix_set(ref_coords, i, j, gsl_matrix_get( svd_coords, it->first, j ) + gsl_matrix_get( ref_coords, i, j));
				}
				++tally;
			}

		}
		// then normalize the the addition
		if (tally != 0)
			for( unsigned int j = 0; j < mColumns; ++j )
		{
			gsl_matrix_set(ref_coords, i, j, gsl_matrix_get( ref_coords, i, j )/((float) tally));
		}
	}

	vector<vector<float> > distances = getDistances(svd_coords, ref_coords);

#if SVD_CLASSIFIER_DEBUG
	cout << "The distance matrix is " << endl;
	for( unsigned int i = 0; i < distances.size(); ++i )
	{
		for( unsigned int j = 0; j < distances[i].size(); ++j )
		{
			cout << distances[i][j] << " ";
		}
		cout << endl;
	}
#endif


	// Finally decide which of the randomly chosen vectors is closest to each of the input vectors
	// and use that as the basis of the grouping
	map< unsigned int, unsigned int> return_map = getMapping(distances);

#if SVD_CLASSIFIER_DEBUG
	cout << "Printing classification map" << endl;
	print_map(return_map);
#endif

	gsl_matrix_free(ref_coords);

	return return_map;
}


map< unsigned int, unsigned int> BoxSVDClassifier::randomSeedCluster(const gsl_matrix* const svd_coords, unsigned int matrix_dims)
{
	// Seed the random number generator
	srand(static_cast<unsigned int>(time(0)));

	vector<unsigned int> random_seed_indices;
	while ( random_seed_indices.size() < mClasses )
	{
		unsigned int random_idx = static_cast<int>(((float)rand()/RAND_MAX)*matrix_dims);
		if ( find( random_seed_indices.begin(), random_seed_indices.end(), random_idx ) == random_seed_indices.end() )
		{
			random_seed_indices.push_back( random_idx );
		}
	}

	// Space to store the reference vectors
	gsl_matrix * ref_coords = gsl_matrix_calloc( mClasses, mColumns );

	// Put the reference vectors into a matrix to make the approach transparent to the reader
	for(unsigned int i = 0; i < random_seed_indices.size(); ++i)
	{
		for( unsigned int j = 0; j < matrix_dims; ++j )
		{
			gsl_matrix_set(ref_coords, i, j, gsl_matrix_get( svd_coords, random_seed_indices[i], j ));
		}
	}

#if SVD_CLASSIFIER_DEBUG
	printMatrix( ref_coords, mClasses, matrix_dims, "Reference matrix in first grouping");
#endif

	// accrue the distance data - this could be done more concisely, but there shouldn't be much cost
	// because the data should be fairl small. By more concisely I mean, the distance data would not need
	// to be stored, it could be determined without storing it in distances.
	vector<vector<float> > distances = getDistances(svd_coords, ref_coords);

#if SVD_CLASSIFIER_DEBUG
	cout << "The distance matrix is " << endl;
	for( unsigned int i = 0; i < distances.size(); ++i )
	{
		for( unsigned int j = 0; j < distances[i].size(); ++j )
		{
			cout << distances[i][j] << " ";
		}
		cout << endl;
	}
#endif


	// Finally decide which of the randomly chosen vectors is closest to each of the input vectors
	// and use that as the basis of the grouping
	map< unsigned int, unsigned int> return_map = getMapping(distances);

#if SVD_CLASSIFIER_DEBUG
	cout << "Printing classification map, randomly seeded" << endl;
	print_map(return_map);
#endif

	gsl_matrix_free(ref_coords);

	return return_map;
}


vector<vector<float> > BoxSVDClassifier::getDistances( const gsl_matrix* const svd_coords, const gsl_matrix* const ref_coords)
{
	// accrue the distance data - this could be done more concisely, but there shouldn't be much cost
	// because the data should be fairl small. By more concisely I mean, the distance data would not need
	// to be stored, it could be determined without storing it in distances.
	vector<vector<float> > distances;
	for (unsigned int i = 0; i < mColumns; ++i )
	{
		vector<float> ith_distances;
		for( unsigned int random_seed_idx = 0; random_seed_idx < mClasses; ++random_seed_idx )
		{
			float distance = 0;
			for (unsigned int j = 0; j < mColumns; ++j )
			{
				float value = (float)( (gsl_matrix_get( ref_coords, random_seed_idx, j) - gsl_matrix_get( svd_coords, i, j)) );
				distance += value * value;
			}
			ith_distances.push_back(distance);
		}
		distances.push_back(ith_distances);
	}

	return distances;
}

map< unsigned int, unsigned int> BoxSVDClassifier::getMapping(const vector<vector<float> >& distances)
{
	// Finally decide which of the randomly chosen vectors is closest to each of the input vectors
	// and use that as the basis of the grouping
	map< unsigned int, unsigned int> return_map;
	unsigned int vector_idx = 0;
	for( vector<vector<float> >::const_iterator it = distances.begin(); it != distances.end(); ++it, ++vector_idx )
	{
		vector<float>::const_iterator mIt = it->begin();
		float min = *mIt;
		unsigned int min_idx = 0;
		for ( unsigned int current_idx = 0; mIt != it->end(); ++mIt, ++current_idx )
		{
			if ( *mIt < min )
			{
				min = *mIt;
				min_idx = current_idx;
			}
		}
		return_map[vector_idx] = min_idx;
	}

	return return_map;
}

map< unsigned int, unsigned int> BoxSVDClassifier::colorMappingByClassSize( const map< unsigned int, unsigned int>& grouping )
{

	vector<unsigned int> current_mappings;
	// Get the extent of the current mappings
	for (map< unsigned int, unsigned int>::const_iterator it = grouping.begin(); it != grouping.end(); ++it )
	{
		if ( find( current_mappings.begin(), current_mappings.end(), it->second ) == current_mappings.end() )
		{
			current_mappings.push_back( it->second );
		}
	}

	if ( current_mappings.size() < 2 )
	{
		cerr << "Error, cannot call colMappingByClassSize when less than 2 classes have been specified, I think you created " << current_mappings.size() << " classes " << endl;
		throw;
	}

	// Record how many data points are in each class.
	map<unsigned int, unsigned int> mappings_tally;
	for( vector<unsigned int>::const_iterator it = current_mappings.begin(); it != current_mappings.end(); ++it )
	{
		// First initialize each total to zero
		mappings_tally[*it] = 0;
	}

	// Now do the actual counting
	for (map< unsigned int, unsigned int>::const_iterator it = grouping.begin(); it != grouping.end(); ++it )
	{
		mappings_tally[it->second] += 1;
	}

	// find the largest tally
	unsigned int current_mapping_idx = 0;
	map< unsigned int, unsigned int> return_map;
	while ( mappings_tally.size() > 0 )
	{
#if SVD_CLASSIFIER_DEBUG
		cout << "Printing mappings_tally" << endl;
		print_map(mappings_tally);
#endif

		map< unsigned int, unsigned int>::iterator it = mappings_tally.begin();
		map< unsigned int, unsigned int>::iterator track_it = mappings_tally.begin();
		unsigned int current_max = it->second;
		unsigned int current_idx = it->first;
		++it;
		for (; it != mappings_tally.end(); ++it )
		{
			if ( it->second > current_max )
			{
				current_max = it->second;
				current_idx = it->first;
				track_it = it;
			}
		}

#if SVD_CLASSIFIER_DEBUG
		cout << "The mapping is " << current_idx << " to " << current_mapping_idx << endl;
#endif
		for (map< unsigned int, unsigned int>::const_iterator group_it = grouping.begin(); group_it != grouping.end(); ++group_it )
		{
			if ( group_it->second == current_idx )
			{
				return_map[group_it->first] = current_mapping_idx;
			}
		}

		mappings_tally.erase( current_idx );

		current_mapping_idx++;
	}


#if SVD_CLASSIFIER_DEBUG
	cout << "Printing adjusted classification map" << endl;
	print_map(return_map);
#endif


	return return_map;
}



vector<float> BoxingTools::get_min_delta_profile(const EMData* const image, int x, int y, int radius)
{
	float peakval = image->get_value_at(x,y);

	vector<float> profile(radius,0); // this sets the vectors size to radius, and the values to 0
	int radius_squared = radius*radius;

	static vector<float> squared_numbers;
	if ( (unsigned int)(radius+1) > squared_numbers.size() ) {
		for(int i = squared_numbers.size(); i <= radius; ++i) {
			squared_numbers.push_back((float)(i*i));
		}
	}

	vector<unsigned int> tally;
	if (mode == SWARM_AVERAGE_RATIO) tally.resize(profile.size(),0);

	for(int k = -radius; k <= radius; ++k) {
		for(int j = -radius; j <= radius; ++j) {
			// Translate coordinates
			int xx = x+j;
			int yy = y+k;

			// Protect against accessing pixels out of bounds
			if ( xx >= image->get_xsize() || xx < 0 ) continue;
			if ( yy >= image->get_ysize() || yy < 0 ) continue;

			// We don't need to pay attention to the origin
			if ( xx == x && yy == y) continue;

			// Protect against vector accesses beyond the boundary
			int square_length = k*k + j*j;
			if (square_length > radius_squared ) continue;

			// The idx is the radius, rounded down. This creates a certain type of pattern that
			// can only really be explained visually...
			int idx = -1;
			// This little loop avoids the use of sqrtf
			for(int i = 1; i < radius; ++i) {
				if ( square_length >= squared_numbers[i] && square_length <= squared_numbers[i+1] ) {
					idx = i;
				}
			}
// 			int idx = (int) sqrtf(k*k + j*j);
			// decrement the idx, because the origin information is redundant
			idx -= 1;

			if ( mode == SWARM_DIFFERENCE ) {
				// Finally, get the drop relative to the origin
				float val = peakval - image->get_value_at(xx,yy);

				// Store it if the drop is smaller than the current value (or if there is no value)
				if ( profile[idx] > val || profile[idx] == 0 ) profile[idx] = val;
			}
			else if (mode == SWARM_RATIO) {
				// Finally, get the drop relative to the origin
				float val =  (peakval - image->get_value_at(xx,yy) ) / peakval;

				// Store it if the drop is smaller than the current value (or if there is no value)
				if ( profile[idx] > val || profile[idx] == 0 ) profile[idx] = val;
			}
			else if (mode == SWARM_AVERAGE_RATIO) {
				profile[idx] += image->get_value_at(xx,yy);
				tally[idx]++;
			}

		}
	}

	if (mode == SWARM_AVERAGE_RATIO) {
		for(unsigned int i = 0; i < profile.size(); ++i) {
			if (tally[i] != 0) {
				profile[i] /= static_cast<float>(tally[i]);
				profile[i] = (peakval - profile[i] ) / peakval;
			}
		}
	}

	return profile;
}

bool BoxingTools::is_local_maximum(const EMData* const image, int x, int y, int radius,EMData* exclusion_map)
{
	float peakval = image->get_value_at(x,y);
	int radius_squared = radius*radius;
	for(int k = -radius; k <= radius; ++k) {
		for(int j = -radius; j <= radius; ++j) {
			// Translate coordinates
			int xx = x+j;
			int yy = y+k;

// 			// Protect against accessing pixel out of bounds
			if ( xx >= image->get_xsize() || xx < 0 ) continue;
			if ( yy >= image->get_ysize() || yy < 0 ) continue;

			// We don't need to pay attention to the origin
			if ( xx == x && yy == y) continue;

			if ((k*k+j*j)>radius_squared) continue;

			if ( image->get_value_at(xx,yy) > peakval)  return false;
		}
	}

	set_radial_non_zero(exclusion_map,x,y,radius);

	return true;

}

vector<IntPoint> BoxingTools::auto_correlation_pick(const EMData* const image, float threshold, int radius, const vector<float>& profile, EMData* const exclusion, const int cradius, int mode)
{
	if (mode < 0 || mode > 2 ) {
		throw InvalidValueException(mode,"Error, the mode can only be 0,1, or 2.");
	}

	if ( radius < 0) {
		throw InvalidValueException(radius,"Radius must be greater than 1");
	}

	if ( cradius < 0) {
		throw InvalidValueException(cradius,"CRadius must be greater than 1");
	}


	int nx = image->get_xsize();
	int ny = image->get_ysize();

	vector<IntPoint> solution;

	int r = radius+1;

	for(int j = r; j < ny-r;++j) {
		for(int k = r; k < nx-r;++k) {

			if (exclusion->get_value_at(k,j) != 0 ) continue;

			if (image->get_value_at(k,j) < threshold) continue;
			if ( mode == 0 ) {
				solution.push_back(IntPoint(k,j));
				set_radial_non_zero(exclusion,k,j,radius);
				continue;
			}

			vector<float> p(r,0);

			if (hi_brid(image,k,j,r,exclusion,p)) {
				if ( mode == 1 ) {
					if (p[cradius] >= profile[cradius]) {
						solution.push_back(IntPoint(k,j));
						set_radial_non_zero(exclusion,k,j,radius);
						continue;
					}
				}
				else /* mode == 2 */{
					bool bad = false;
					for (int ii = 0; ii <= cradius; ++ii) {
						if (p[ii] < profile[ii]) {
							bad = true;
							break;
						}
					}
					if (bad) continue;
					solution.push_back(IntPoint(k,j));
					set_radial_non_zero(exclusion,k,j,radius);
				}


			}
		}
	}
	return solution;
}


bool BoxingTools::hi_brid(const EMData* const image, int x, int y, int radius,EMData* const exclusion_map, vector<float>& profile)
{

	float peakval = image->get_value_at(x,y);

	int radius_squared = radius*radius;

	static vector<float> squared_numbers;
	if ( (unsigned int)(radius+1) > squared_numbers.size() ) {
		for(int i = squared_numbers.size(); i <= radius; ++i) {
			squared_numbers.push_back((float)(i*i));
		}
	}

	vector<unsigned int> tally;
	if (mode == SWARM_AVERAGE_RATIO) tally.resize(profile.size(),0);

	for(int k = -radius; k <= radius; ++k) {
		for(int j = -radius; j <= radius; ++j) {
			// Translate coordinates
			int xx = x+j;
			int yy = y+k;

			// Protect against accessing pixels out of bounds
			if ( xx >= image->get_xsize() || xx < 0 ) continue;
			if ( yy >= image->get_ysize() || yy < 0 ) continue;

			// We don't need to pay attention to the origin
			if ( xx == x && yy == y) continue;

			// Protect against vector accesses beyond the boundary
			int square_length = k*k + j*j;
			if (square_length > radius_squared ) continue;

			// It's not a local maximum!
			if ( image->get_value_at(xx,yy) > peakval)  return false;

			// The idx is the radius, rounded down. This creates a certain type of pattern that
			// can only really be explained visually...
			int idx = -1;
			// This little loop avoids the use of sqrtf
			for(int i = 1; i < radius; ++i) {
				if ( square_length >= squared_numbers[i] && square_length <= squared_numbers[i+1] ) {
					idx = i;
				}
			}
// 			int idx = (int) sqrtf(k*k + j*j);
			// decrement the idx, because the origin information is redundant
			idx -= 1;

			if (mode == SWARM_DIFFERENCE) {
				// Finally, get the drop relative to the origin
				float val = peakval - image->get_value_at(xx,yy);

				// Store it if the drop is smaller than the current value (or if there is no value)
				if ( profile[idx] > val || profile[idx] == 0 ) profile[idx] = val;
			}
			else if (mode == SWARM_RATIO) {
				// Finally, get the drop relative to the origin
				float val =  (peakval - image->get_value_at(xx,yy) ) / peakval;

				// Store it if the drop is smaller than the current value (or if there is no value)
				if ( profile[idx] > val || profile[idx] == 0 ) profile[idx] = val;
			}
			else if (mode == SWARM_AVERAGE_RATIO) {
				profile[idx] += image->get_value_at(xx,yy);
				tally[idx]++;
			}

		}
	}

	if (mode == SWARM_AVERAGE_RATIO) {
		for(unsigned int i = 0; i < profile.size(); ++i) {
			if (tally[i] != 0) {
				profile[i] /= static_cast<float>(tally[i]);
				profile[i] = (peakval - profile[i] ) / peakval;
			}
		}
	}

	set_radial_non_zero(exclusion_map,x,y,radius);

	return true;
}


void BoxingTools::set_radial_non_zero(EMData* const exclusion, int x, int y, int radius)
{
	int radius_squared = radius*radius;
	for(int k = -radius; k <= radius; ++k) {
		for(int j = -radius; j <= radius; ++j) {
			// Translate coordinates
			int xx = x+j;
			int yy = y+k;

			if ((k*k+j*j)>radius_squared) continue;
			// Protect against accessing pixel out of bounds
			if ( xx >= exclusion->get_xsize() || xx < 0 ) continue;
			if ( yy >= exclusion->get_ysize() || yy < 0 ) continue;

			exclusion->set_value_at(xx,yy,1);
		}
	}
}

IntPoint BoxingTools::find_radial_max(const EMData* const map, int x, int y, int radius)
{
	float currentmax = map->get_value_at(x,y);

	IntPoint soln(x,y);

	int radius_squared = radius*radius;
	for(int k = -radius; k <= radius; ++k) {
		for(int j = -radius; j <= radius; ++j) {
			// Translate coordinates
			int xx = x+j;
			int yy = y+k;

			// Protect against accessing pixels out of bounds
			if ( xx >= map->get_xsize() || xx < 0 ) continue;
			if ( yy >= map->get_ysize() || yy < 0 ) continue;

			// Protect against vector accesses beyond the boundary
			int square_length = k*k + j*j;
			if (square_length > radius_squared ) continue;

			float val = map->get_value_at(xx,yy);

			if (val > currentmax) {
				currentmax = val;
				soln[0] = xx;
				soln[1] = yy;
			}
		}
	}

	return soln;
}


void BoxingTools::set_region( EMData* const image, const EMData* const mask, const int x, const int y, const float& val) {

	// Works only in 2D
	int inx = image->get_xsize();
	int iny = image->get_ysize();
	int mnx = mask->get_xsize();
	int mny = mask->get_ysize();

	int startx = x-mnx/2;
	int endx =startx + mnx;
	int xoffset = 0;
	if (startx < 0) {
		xoffset = abs(startx);
		startx = 0;
	}
	if (endx > inx) endx = inx;

	int starty = y-mny/2;
	int endy =starty + mny;
	int yoffset = 0;
	if (starty < 0) {
		yoffset = abs(starty);
		starty = 0;
	}
	if (endy > iny) endy = iny;


	for (int j = starty; j < endy; ++j ) {
		for (int i = startx; i < endx; ++i) {
			if (mask->get_value_at(xoffset+i-startx,yoffset+j-starty) != 0 ) {
				image->set_value_at(i,j,val);
			}
		}
	}
}

map<unsigned int, unsigned int> BoxingTools::classify(const vector<vector<float> >& data, const unsigned int& classes)
{
	BoxSVDClassifier classifier(data, classes);
	map< unsigned int, unsigned int> mapping = classifier.go();

	mapping = BoxSVDClassifier::colorMappingByClassSize( mapping );

	return mapping;

}


Vec3f BoxingTools::get_color( const unsigned int index )
{
	if ( colors.size() == 0 ) {
		colors.push_back(Vec3f(0,0,0));
		colors.push_back(Vec3f(0,0,1));
		colors.push_back(Vec3f(0,1,0));
		colors.push_back(Vec3f(1,0,0));
		colors.push_back(Vec3f(1,1,0));
		colors.push_back(Vec3f(1,0,1));
		colors.push_back(Vec3f(0,1,1));
		colors.push_back(Vec3f(1,1,1));
	}
	if ( index >= colors.size() )
	{
		while ( colors.size() <= index )
		{
			bool found = false;
			while ( !found )
			{
				unsigned int random_idx = rand() % colors.size();
				unsigned int random_idx2 = rand() % colors.size();
				while ( random_idx2 == random_idx )
				{
					random_idx2 = rand() % colors.size();
				}

				Vec3f result = (colors[random_idx] + colors[random_idx2])/2.0;
				if ( find( colors.begin(), colors.end(), result ) == colors.end() )
				{
					colors.push_back( result );
					found = true;
				}
			}
		}
	}
	return colors[index];
}

