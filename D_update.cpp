#define _USE_MATH_DEFINES
#include <cmath>
#include "fdtd3d.h"

void E_update( 
    double ***nDr, double ***nDth, double ***nDphi,
    double ***oDr, double ***oDth, double ***oDphi,
    double ***Hr, double ***Hth, double ***Hphi )
    {
        double CDr1, CDr2;
        double CDth1, CDth2;
        double CDph1, CDph2;
        
        double ri1, ri2, ri3;
        double sin_th1, sin_th2, sin_th3;

        for( int i = 0; i < Nr; i++ ){
            for( int j = L + 1; j < Ntheta - L; j++ ){
                for( int k = L + 1; k < Nphi - L; k++ ){
                    CDr1 = Dt/dist(i+0.5)/std::sin(th(j))/delta_theta;
                    CDr2 = Dt/dist(i)/std::sin(th(j+0.5))/delta_phi;

                    nDr[i][j][k] = oDr[i][j][k] + CDr1*( std::sin(th(j+0.5))*Hphi[i][j][k] - std::sin(th(j-0.5))*Hphi[i][j-1][k] )
                                - CDr2*( Hth[i][j][k] - Hth[i][j][k-1] );
                }
            }
        }

        for( int i = 1; i < Nr; i++ ){
            for( int j = L; j < Ntheta - L; j++ ){
                for( int k = L + 1; k < Nphi - L; k++ ){
                    CDth1 = Dt/dist(i)/std::sin(th(j+0.5))/delta_phi;
                    CDth2 = Dt/dist(i)/delta_r;

                    nDth[i][j][k] = oDth[i][j][k] + CDth1*( Hr[i][j][k] - Hr[i][j][k-1] )
                                - CDth2*( dist(i+0.5)*Hth[i][j][k] - dist(i-0.5)*Hth[i-1][j][k] );
                }
            }
        }

        for( int i = 1; i < Nr; i++ ){
            for( int j = L + 1; j < Ntheta - L; j++ ){
                for( int k = L; k < Nphi - L; k++ ){
                    CDph1 = Dt/dist(i)/delta_r;
                    CDph2 = Dt/dist(i)/delta_theta;

                    nDphi[i][j][k] = oDphi[i][j][k] + CDph1*( dist(i+0.5)*Hth[i][j][k] - dist(i-0.5)*Hth[i-1][j][k] )
                                - CDph2*( Hr[i][j][k] - Hr[i][j-1][k] );
                }
            }
        }

    }