#define _USE_MATH_DEFINES
#include <cmath>
#include "fdtd3d.h"

void H_update_pml(
    double ***Er, double ***Eth, double ***Eph,
    double ***Hr, double ***Hth, double ***Hph,
    double ****Hr_th1, double ****Hr_th2, double ****Hr_ph,
    double ****Hth_ph, double ****Hth_r,
    double ****Hph_r, double ****Hph_th,
    double *sigma_th_h, double *sigma_ph_h,
    pml *idx_r, pml *idx_th, pml *idx_ph
){
    for( int area = 0; area < 4; area++ ){
        for( int i = 0; i < Nr + 1; i++ ){
            for( int j = idx_r[area].j1(); j <= idx_r[area].j2(); j++ ){
                int j_area = j - idx_r[area].j1();
                for( int k = idx_r[area].k1(); k <= idx_r[area].k2(); k++ ){
                    int k_area = k - idx_r[area].k1();

                    Hr_th1[area][i][j_area][k_area] = C_1(sigma_th_h[j])*Hr_th1[area][i][j_area][k_area]
                        - C_2(dist(i), sigma_th_h[j])/MU0*(Eph[i][j+1][k] - Eph[i][j][k]);

                    Hr_th2[area][i][j_area][k_area] = Hr_th2[area][i][j_area][k_area]
                        - C_3(dist(i), th(j+0.5))/MU0*(Eph[i][j+1][k] + Eph[i][j][k]);
                    
                    Hr_ph[area][i][j_area][k_area] = C_1(sigma_ph_h[k])*Hr_ph[area][i][j_area][k_area]
                        + C_4(dist(i), th(j+0.5), sigma_ph_h[k])/MU0*(Eth[i][j][k+1] - Eth[i][j][k]);

                    Hr[i][j][k] = Hr_th1[area][i][j_area][k_area]
                            + Hr_th2[area][i][j_area][k_area] + Hr_ph[area][i][j_area][k_area];
                }
            }
        }
    }

    for(int area = 0; area < 4; area++ ){
        for( int i = 0; i < Nr; i++ ){
            for( int j = idx_th[area].j1(); j <= idx_th[area].j2(); j++ ){
                int j_area = j - idx_th[area].j1();
                for( int k = idx_th[area].k1(); k <= idx_th[area].k2(); k++ ){
                    int k_area = k - idx_th[area].k1();

                    Hth_ph[area][i][j_area][k_area] = C_1( sigma_ph_h[k] )*Hth_ph[area][i][j_area][k_area]
                                - C_4(dist(i+0.5), th(j), sigma_ph_h[k])/MU0*(Er[i][j][k+1] - Er[i][j][k]);
                    Hth_r[area][i][j_area][k_area] = Hth_r[area][i][j_area][k_area]
                                + C_5(dist(i+0.5))/MU0*(dist(i+1.0)*Eph[i+1][j][k] - dist(i)*Eph[i][j][k]);
                    
                    Hth[i][j][k] = Hth_ph[area][i][j_area][k_area] + Hth_r[area][i][j_area][k_area];

                }
            }
        }
    }

    for( int area = 0; area < 4; area++ ){
        for( int i = 0; i < Nr; i++ ){
            for( int j = idx_ph[area].j1(); j <= idx_ph[area].j2(); j++ ){
                int j_area = j - idx_ph[area].j1();
                for( int k = idx_ph[area].k1(); k <= idx_ph[area].k2(); k++ ){
                    int k_area = k - idx_ph[area].k1();

                    Hph_r[area][i][j_area][k_area] = Hph_r[area][i][j_area][k_area]
                                - C_5(dist(i+0.5))/MU0*(dist(i+1.0)*Eth[i+1][j][k] - dist(i)*Eth[i][j][k]);
                    
                    Hph_th[area][i][j_area][k_area] = C_1(sigma_th_h[j])*Hph_th[area][i][j_area][k_area]
                                + C_6(dist(i+0.5), sigma_th_h[j])/MU0*(Er[i][j+1][k] - Er[i][j][k]);
                    
                    Hph[i][j][k] = Hph_r[area][i][j_area][k_area] + Hph_th[area][i][j_area][k_area];
                }
            }
        }
    }
}