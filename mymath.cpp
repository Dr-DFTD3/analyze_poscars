#include "mymath.h"


double** 
invert_matrix(double** H)
{
    double** Hinv = new double*[3];
    for (int i=0;i<3;i++) Hinv[i] = new double[3];

    double detm =  H[0][0]*H[1][1]*H[2][2] - H[0][0]*H[2][1]*H[1][2] \
        - H[1][0]*H[0][1]*H[2][2] + H[1][0]*H[2][1]*H[0][2] \
        + H[2][0]*H[0][1]*H[1][2] - H[2][0]*H[1][1]*H[0][2];

    Hinv[0][0] = (H[2][2]*H[3][3] - H[2][3]*H[3][2])/detm;
    Hinv[1][0] = (H[2][3]*H[3][1] - H[2][1]*H[3][3])/detm;
    Hinv[2][0] = (H[2][1]*H[3][2] - H[2][2]*H[3][1])/detm;
    
    Hinv[0][1] = (H[3][2]*H[1][3] - H[3][3]*H[1][2])/detm;
    Hinv[1][1] = (H[3][3]*H[1][1] - H[3][1]*H[1][3])/detm;
    Hinv[2][1] = (H[3][1]*H[1][2] - H[3][2]*H[1][1])/detm;
    
    Hinv[0][2] = (H[1][2]*H[2][3] - H[1][3]*H[2][2])/detm;
    Hinv[1][2] = (H[1][3]*H[2][1] - H[1][1]*H[2][3])/detm;
    Hinv[2][2] = (H[1][1]*H[2][2] - H[1][2]*H[2][1])/detm;

    return Hinv;
}