#define _USE_MATH_DEFINES
#include <cmath>
#include "fdtd3d.h"

void D_update_pml(
    double ***Dr, double ***Dth, double ***Dph,
    double ***Hr, double ***Hth, double ***Hph,
    double ****Dr_th1, double ****Dr_th2, double ****Dr_ph,
    double ****Dth_ph, double ****Dth_r,
    double ****Dph_r, double ****Dph_th,
    double *sigma_th, double *sigma_ph,
    pml *idx_r, pml *idx_th, pml *idx_ph)
    {
        for( int area = 0; area < 4; area++ ){
            for(int i = 0; i < Nr; i++ ){
                for( int j = idx_r[area].j1(); j <=idx_r[area].j2(); j++ ){
                    int j_area = j - idx_r[area].j1();
                    for( int k = idx_ph[area].k1(); k <= idx_ph[area].k2(); k++ ){
                        int k_area = k - idx_r[area].k1();
                        Dr_th1[area][i][j_area][k_area] = C_1(sigma_th[j])*Dr_th1[area][i][j_area][k_area]
                                                    + C_2(dist(i + 0.5), sigma_th[j])*(Hph[i][j][k] - Hph[i][j-1][k]);
                        Dr_th2[area][i][j_area][k_area] = Dr_th2[area][i][j_area][k_area]
                                                    + C_3(dist(i+0.5), th(j))*(Hph[i][j][k] + Hph[i][j-1][k]);
                        Dr_ph[area][i][j_area][k_area] = C_1(sigma_ph[k])*Dr_ph[area][i][j_area][k_area]
                                                    - C_4(dist(i+0.5), th(j), sigma_ph[k])*(Hth[i][j][k] - Hth[i][j][k-1]);

                        Dr[i][j][k] = Dr_th1[area][i][j_area][k_area]
                                + Dr_th2[area][i][j_area][k_area] + Dr_ph[area][i][j_area][k_area];
                    }
                }
            }
        }

        for( int area = 0; area < 4; area++ ){
            for( int i = 1; i < Nr; i++ ){
                for( int j = idx_th[area].j1(); j <= idx_th[area].j2(); j++ ){
                    int j_area = j - idx_th[area].j1();
                    for( int k = idx_ph[area].k1(); k <= idx_ph[area].k2(); k++ ){
                        int k_area = k - idx_ph[area].k1();

                        Dth_ph[area][i][j_area][k_area] = C_1(sigma_ph[k])*Dth_ph[area][i][j_area][k_area]
                                                    + C_4(dist(i), th(j+0.5), sigma_ph[k])*(Hr[i][j][k] - Hr[i][j][k-1]);
                        Dth_r[area][i][j_area][k_area] = Dth_r[area][i][j_area][k_area]
                                                    - C_5(dist(i))*(dist(i+0.5)*Hph[i][j][k] - dist(i-0.5)*Hph[i-1][j][k]);
                        
                        Dth[i][j][k] = Dth_ph[area][i][j_area][k_area] + Dth_r[area][i][j_area][k_area];

                    }
                }
            }
        }

        for( int area = 0; area < 4; area++ ){
            for( int i = 1; i < Nr; i++ ){
                for( int j = idx_ph[area].j1(); j <= idx_ph[area].j2(); j++ ){
                    int j_area = j - idx_ph[area].j1();
                    for( int k = idx_ph[area].k1(); k <= idx_ph[area].k2(); k++ ){
                        int k_area = idx_ph[area].k1();
                        Dph_r[area][i][j_area][k_area] = Dph_r[area][i][j_area][k_area]
                                                    + C_5(dist(i))*(dist(i+0.5)*Hth[i][j][k] - dist(i-0.5)*Hth[i][j-1][k]);
                        Dph_th[area][i][j_area][k_area] = C_1(sigma_th[j])*Dph_th[area][i][j_area][k_area]
                                                    - C_6(dist(i), sigma_th[j])*( Hr[i][j][k] - Hr[i][j-1][k] );
                        
                        Dph[i][j][k] = Dph_r[area][i][j_area][k_area] + Dph_th[area][i][j_area][k_area];
                    }
                }
            }
        }

    }