#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <chrono>
#include <Eigen/Core>
#include "fdtd3d.h"

const int Nr{100};
const int Ntheta{100};
const int Nphi{1000};

constexpr double R_r{100.0e3};

const double delta_r{ R_r/(double)Nr };
const double delta_theta{ 1.0e3/(double)R0 };
const double delta_phi{ 1.0e3/(double)R0 };
const double Dt { double( 0.99/C0/std::sqrt(1.0/delta_r/delta_r
 + 1.0/R0/R0/delta_theta/delta_theta
 + 1.0/R0/R0/std::sin(THETA0)/std::sin(THETA0)/delta_phi/delta_phi) ) };
const double inv_Dt{ 1.0/Dt };
const double sigma_t{ 7.0*Dt };
const double t0{ 6.0*sigma_t };

// center point //
const int i_0{ Nr/2 };
const int j_0{ Ntheta/2 };
const int k_0{ Nphi/2 };

// Source point //
const int i_s{1};
const int j_s{50};
const int k_s{100};

// Receive Point //
const int i_r{1};
const int j_r{50};
const int k_r{ Nphi - 50 };

// PML info //
const int L{10};
const double M{3.5};
const double R{1.0e-6};

const double sigma_th_max{ -(M + 1.0)*C0*std::log(R)/2.0/double(L)/delta_theta/R0 };
const double sigma_phi_max{ -(M + 1.0)*C0*std::log(R)/2.0/double(L)/delta_phi/R0 }

int main( void )
{
    int time_step = 1700;
    double t;
    double J;
    
    double ***Hr, ***Htheta, ***Hphi;
    Hr = memory_allocate( Nr+1, Ntheta, Nphi, 0.0 );
    Htheta = memory_allocate( Nr, Ntheta+1, Nphi, 0.0 );
    Hphi = memory_allocate( Nr, Ntheta, Nphi+1, 0.0 );

    double ***nEr, ***nEtheta, ***nEphi;
    nEr = memory_allocate3d( Nr, Ntheta+1, Nphi+1, 0.0 );
    nEtheta = memory_allocate3d( Nr+1, Ntheta, Nphi+1, 0.0 );
    nEphi = memory_allocate3d( Nr+1, Ntheta+1, Nphi, 0.0 );

    double ***oEr, ***oEtheta, ***oEphi;
    oEr = memory_allocate3d( Nr, Ntheta+1, Nphi+1, 0.0 );
    oEtheta = memory_allocate3d( Nr+1, Ntheta, Nphi+1, 0.0 );
    oEphi = memory_allocate3d( Nr+1, Ntheta+1, Nphi, 0.0 );

    double ***Dr, ***Dtheta, ***Dphi;
    Dr = memory_allocate3d( Nr, Ntheta+1, Nphi+1, 0.0 );
    Dtheta = memory_allocate3d( Nr+1, Ntheta, Nphi+1, 0.0 );
    Dphi = memory_allocate3d( Nr+1, Ntheta+1, Nphi, 0.0 );

    double *sigma_theta, *sigma_phi, *sigma_theta_h, *sigma_phi_h;
    sigma_theta = new double[Ntheta + 1];
    sigma_phi = new double[Nphi + 1];
    sigma_theta_h = new doubke[Ntheta + 1];
    sigma_phi_h = new double[Nphi + 1];

}
