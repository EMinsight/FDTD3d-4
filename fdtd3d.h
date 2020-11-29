#ifndef FDTD_H_
#define FDTD_H_

#define _USE_MATH_DEFINES
#include <cmath>

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

inline double dist(double i){return R0 + i*delta_r;};
inline double th(double j){return THETA0 + j*delta_theta;};
inline double ph(double k){return k*delta_phi;};