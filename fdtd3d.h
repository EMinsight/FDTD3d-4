#ifndef FDTD_H_
#define FDTD_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include "pml.h"

#define C0 (3.0e8)
#define MU0 (4.0*M_PI*1.0e-7)
#define EPS0 (1.0/MU0/C0/C0)
#define R0 (6370e3)
#define THETA0 (M_PI*0.5 - std::atan(50e3/R0))
#define E_Q (1.6e-19)
#define E_M (9.11e-31)

extern const int Nr;
extern const int Ntheta;
extern const int Nphi;

extern const double delta_r;
extern const double delta_theta;
extern const double delta_phi;
extern const double Dt;
extern const double inv_Dt;

extern const int L;
extern const double M;
extern const double R;

extern const double sigma_th_max;
extern const double sigma_phi_max;

double*** memory_allocate3d( int, int, int, double );

void sigma_calc( double *sigma_th, double *sigma_ph, double *sigma_th_h, double *sigma_ph_h );
void PML_idx_initialize(
                        pml* idx_Dr, pml* idx_Dth, pml*idx_Dphi,
                        pml* idx_Hr, pml* idx_Hth, pml* idx_Hphi);

void PML_field_initialize(
    double**** Dr_th1, double**** Dr_th2, double**** Dr_ph,
    double**** Dth_ph, double**** Dth_r,
    double**** Dph_r, double**** Dph_th,
    double**** Hr_th1, double**** Hr_th2, double**** Hr_ph,
    double**** Hth_ph, double**** Hth_r,
    double**** Hph_r, double**** Hph_th
);

void D_update( double ***newDr, double ***newDth, double ***newDph,
                double ***oldDr, double ***oldDth, double ***oldDph,
                double ***Hr, double ***Hth, double ***Hph );

void D_update_pml(
    double ***Dr, double ***Dth, double ***Dph, double ***Hr, double ***Hth, double ***Hph,
    double ****Dr_th1, double ****Dr_th2, double ****Dr_ph,
    double ****Dth_ph, double ****Dth_r,
    double ****Dph_r, double ****Dph_th,
    double *sigma_th, double *sigma_ph,
    pml *idx_Dr, pml *idx_Dth, pml *idx_Dph
);
void E_update(
    double ***nEr, double ***nEth, double ***nEph, double ***oEr, double ***oEth, double ***oEph,
    double ***nDt, double ***nDth, double ***nDph, double ***oDr, double ***oDth, double ***oDph
);

void H_update(
    double ***nEr, double ***nEth, double ***nEph,
    double ***Hr, double ***Hth, double ***Hph
);

void H_update_pml(
    double ***nEr, double ***nEth, double ***nEph, double ***Hr, double ***Hth, double ***Hph,
    double ****Hr_th1, double ****Hr_th2, double ****Hr_ph,
    double ****Hth_ph, double ****Hth_r,
    double ****Hph_r, double ****Hph_th,
    double *sigma_th_h, double *sigma_ph_h, pml *idx_Hr, pml *idx_Hth, pml *idx_Hph
);

inline double dist(double i){return R0 + i*delta_r;};
inline double th(double j){return THETA0 + j*delta_theta;};
inline double ph(double k){return k*delta_phi;};

inline double C_1(double sig){return ((inv_Dt - sig/2.0)/(inv_Dt + sig/2.0));};
inline double C_2(double r, double sig){return 1.0/r/delta_theta/(inv_Dt + sig/2.0);};
inline double C_3(double r, double theta){return Dt*std::cos(theta)/std::sin(theta)/2.0/r;};
inline double C_4(double r, double theta, double sig){return 1.0/r/std::sin(theta)/delta_phi/(inv_Dt + sig/2.0);};
inline double C_5(double r){return Dt/r/delta_r;};
inline double C_6(double r, double sig){return 1.0/(inv_Dt + sig/2.0)/r/delta_theta;};