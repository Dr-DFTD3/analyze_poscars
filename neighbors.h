#ifndef _NEIGHS_H_
#define _NEIGHS_H_
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <fstream>
#include <iostream>
#include <iomanip>

#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "common.h"
#include "mymath.h"


/*------------------------------------------------------------
    Construct a nearest neighbor list for each atom in the 
    system. Use an n^2 recursive search for nearest neighbors 
    within a fixed cutoff RADIUS, or use covalent radii for 
    pair of atoms. Returns nothing, update the neighborlist 
    of each atom
--------------------------------------------------------------*/
void 
build_nn_list( std::vector<point_t>& atoms,
			   double**              H,
			   param_t*              param );


/*------------------------------------------------------------
    Cycle through each atoms' nearest neighbor list to 
    identify clusters of atoms, i.e. molecules.
    Assume each atom is in a cluster by itself, recursivley 
    update which cluster each atom  belongs to. Return 
    an array containing the member ids of each atom belonging
    to a cluster.
--------------------------------------------------------------*/
std::vector<std::vector<int> >
find_clusters( char*                     stdName, 
			   std::vector<point_t>&     atoms, 
			   std::vector<std::string>& species,
			   param_t*                  param,
               std::vector<int>&         clusterTags );


#endif