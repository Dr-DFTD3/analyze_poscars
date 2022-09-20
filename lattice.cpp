#include "lattice.h"
#include "mymath.h"


/*------------------------------------------------------------
    Convert fractional [0,1] coordinates to cartesian, in Å
--------------------------------------------------------------*/
void 
cart_coords( std::vector<point_t>& atoms, 
             double** H)
{
    double xs,ys,zs;

    for (int n=0;n<atoms.size();n++) 
    {
        xs = atoms[n].xs;
        ys = atoms[n].ys;
        zs = atoms[n].zs;

        atoms[n].x = xs*H[0][0] + ys*H[1][0] + zs*H[2][0];
        atoms[n].y = xs*H[0][1] + ys*H[1][1] + zs*H[2][1];
        atoms[n].z = xs*H[0][2] + ys*H[1][2] + zs*H[2][2];
    }
}

/*------------------------------------------------------------
    Convert cartesian coordinates to fractional [0,1] coordinates
--------------------------------------------------------------*/
void 
reduce_coords( std::vector<point_t>& atoms, 
               double** H)
{
    double** Hinv = invert_matrix( H );

    double x,y,z;

    for (int n=0;n<atoms.size();n++) 
    {
        x = atoms[n].x;
        y = atoms[n].y;
        z = atoms[n].z;

        atoms[n].xs = x*Hinv[0][0] + y*Hinv[1][0] + z*Hinv[2][0];
        atoms[n].ys = x*Hinv[0][1] + y*Hinv[1][1] + z*Hinv[2][1];
        atoms[n].zs = x*Hinv[0][2] + y*Hinv[1][2] + z*Hinv[2][2];
    }
}
