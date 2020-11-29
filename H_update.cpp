#define _USE_MATH_DEFINES
#include <cmath>
#include "fdtd3d.h"

void H_update(
    double ***Er, double ***Eth, double ***Eph,
    double ***Hr, double ***Hth, double ***Hph
){
    double CHr1, CHr2;
    double CHth1, CHth2;
    double CHph1, CHph2;

    for( int i = 1; i < Nr + 1; i++ ){
        for( int j = L; j < Ntheta - L; j++ ){
            for( int k = L; k < Nphi - L; k++ ){
                CHr1 = Dt/MU0/dist(i)/std::sin(th(j+0.5))/delta_theta;
                CHr1 = Dt/MU0/dist(i)/std::sin(th(j+0.5))/delta_phi;

                Hr[i][j][k] = Hr[i][j][k] - CHr1*(std::sin(th(j+0.5))*Eph[i][j+1][k] - std::sin(th(j))*Eph[i][j][k])
                    + CHr2*(Eth[i][j][k+1] - Eth[i][j][k]);
            }
        }
    }

    for( int i = 0; i < Nr; i++ ){
        for( int j = L + 1; j < Ntheta - L; j++ ){
            for( int k = L; k < Nphi - L; k++ ){
                CHth1 = Dt/MU0/dist(i+0.5)/std::sin(th(j))/delta_phi;
                CHth2 = Dt/MU0/dist(i+0.5)/delta_r;

                Hth[i][j][k] = Hth[i][j][k] - CHth1*(Er[i][j][k+1] - Er[i][j][k])
                        + CHth2*(dist(i+1.0)*Eph[i+1][j][k] - dist(i)*Eph[i][j][k]);

            }
        }
    }

    for( int i = 0; i < Nr; i++ ){
        for( int j = L; j < Ntheta - L; j++ ){
            for( int k = L + 1; k < Nphi - L; k++ ){
                CHph1 = Dt/MU0/dist(i+0.5)/delta_r;
                CHph2 = Dt/MU0/dist(i+0.5)/delta_theta;

                Hph[i][j][k] = Hph[i][j][k] - CHph1*(dist(i+1.0)*Eth[i+1][j][k] - dist(i)*Eth[i][j][k])
                        + CHph2*(Er[i][j+1][k] - Er[i][j][k]);
            }
        }
    }

}