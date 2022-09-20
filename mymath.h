#ifndef _MYMATH_H_
#define _MYMATH_H_
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


template <class Type> 
const Type& 
max (const Type& a, const Type& b) {
    return (a<b)?b:a;     // or: return comp(a,b)?b:a; for version (2)
}

template <class Type> 
const Type& 
min (const Type& a, const Type& b) 
{
    return !(b<a)?a:b;     // or: return !comp(b,a)?a:b; for version (2)
}

double** 
invert_matrix(double** H);


#endif
