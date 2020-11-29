#define _USE_MATH_DEFINES
#include <cmath>
#include "fdtd3d.h"

void E_update(
    double ***nEr, double ***nEth, double ***nEph,
    double ***oEr, double ***oEth, double ***oEph,
    double ***nDr, double ***nDth, double ***nDph,
    double ***oDr, double ***oDth, double ***oDph
){
    for( int i = 0; i < Nr; i++ ){
        for( int j = 1; j < Ntheta; j++ ){
            for( int k = 1; k < Nphi; k++ ){
                nDr[i][j][k] = oDr[i][j][k] +
                    (nDr[i][j][k] - oDr[i][j][k])/EPS0;
            }
        }
    }

    for( int i = 1; i < Nr; i++ ){
        for( int j = 0; j < Ntheta; j++ ){
            for( int k = 1; k < Nphi; k++ ){
                nEth[i][j][k] = oEth[i][j][k] + 
                    (nDth[i][j][k] - oDth[i][j][k])/EPS0;
            }
        }
    }

    for( int i = 1; i < Nr; i++ ){
        for( int j = 1; j < Ntheta; j++ ){
            for( int k = 0; k < Nphi; k++ ){
                nEph[i][j][k] = oEph[i][j][k] + 
                        (nDph[i][j][k] - oDph[i][j][k])/EPS0;
            }
        }
    }

}