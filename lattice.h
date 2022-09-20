#ifndef _LATTICE_H_
#define _LATTICE_H_
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

/*------------------------------------------------------------
    Convert cartesian coordinates to fractional [0,1] coordinates
--------------------------------------------------------------*/
void 
reduce_coords( std::vector<point_t>& atoms, 
			   double**              H);

/*------------------------------------------------------------
    Convert fractional [0,1] coordinates to cartesian, in Å
--------------------------------------------------------------*/
void 
cart_coords( std::vector<point_t>& atoms, 
	         double**              H);

#endif
