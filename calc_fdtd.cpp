#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
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
const double sigma_ph_max{ -(M + 1.0)*C0*std::log(R)/2.0/double(L)/delta_phi/R0 };

int calc_fdtd( void )
{
    int time_step = 1700;
    double t = 0.0;
    double J = 0.0;
    
    double ***Hr, ***Hth, ***Hph;
    Hr = memory_allocate( Nr+1, Ntheta, Nphi, 0.0 );
    Hth = memory_allocate( Nr, Ntheta+1, Nphi, 0.0 );
    Hph = memory_allocate( Nr, Ntheta, Nphi+1, 0.0 );

    double ***nEr, ***nEth, ***nEph;
    nEr = memory_allocate3d( Nr, Ntheta+1, Nphi+1, 0.0 );
    nEth = memory_allocate3d( Nr+1, Ntheta, Nphi+1, 0.0 );
    nEph = memory_allocate3d( Nr+1, Ntheta+1, Nphi, 0.0 );

    double ***oEr, ***oEth, ***oEph;
    oEr = memory_allocate3d( Nr, Ntheta+1, Nphi+1, 0.0 );
    oEth = memory_allocate3d( Nr+1, Ntheta, Nphi+1, 0.0 );
    oEph = memory_allocate3d( Nr+1, Ntheta+1, Nphi, 0.0 );

    double ***nDr, ***nDth, ***nDph;
    nDr = memory_allocate3d( Nr, Ntheta+1, Nphi+1, 0.0 );
    nDth = memory_allocate3d( Nr+1, Ntheta, Nphi+1, 0.0 );
    nDph = memory_allocate3d( Nr+1, Ntheta+1, Nphi, 0.0 );

    double ***oDr, ***oDth, ***oDph;
    oDr = memory_allocate3d( Nr, Ntheta+1, Nphi+1, 0.0 );
    oDth = memory_allocate3d( Nr+1, Ntheta, Nphi+1, 0.0 );
    oDph = memory_allocate3d( Nr+1, Ntheta+1, Nphi, 0.0 );

    double *sigma_th, *sigma_ph, *sigma_th_h, *sigma_ph_h;
    sigma_th = new double[Ntheta + 1];
    sigma_ph = new double[Nphi + 1];
    sigma_th_h = new doubke[Ntheta + 1];
    sigma_ph_h = new double[Nphi + 1];

    double ****Dr_th1, ****Dr_th2, ****Dr_ph;
    double ****Dth_ph, ****Dth_r;
    double ****Dph_r, ****Dph_th;

    Dr_th1 = new double***[4];
    Dr_th2 = new double***[4];
    Dr_ph = new double***[4];
    Dth_ph = new double***[4];
    Dth_r = new double***[4];
    Dph_r = new double***[4];
    Dph_th = new double***[4];

    double ****Hr_th1, ****Hr_th2, ****Hr_ph;
    double ****Hth_ph, ****Hth_r;
    double ****Hph_r, ****Hph_th;

    Hr_th1 = new double***[4];
    Hr_th2 = new double***[4];
    Hr_ph = new double***[4];
    Hth_ph = new double***[4];
    Hth_r = new double***[4];
    Hph_r = new double***[4];
    Hph_th = new double***[4];

    PML_field_initialize(
        Dr_th1, Dr_th2, Dr_ph,
        Dth_ph, Dph_r,
        Dph_r, Dph_th,
        Hr_th1, Hr_th2, Hr_ph,
        Hth_ph, Hph_r,
        Hph_r, Hph_th
    );

    /*double ***Dr_th1_A1, ***Dr_th2_A1, ***Dr_ph_A1;
    double ***Dr_th1_A2, ***Dr_th2_A2, ***Dr_ph_A2;
    double ***Dr_th1_A3, ***Dr_th2_A3, ***Dr_ph_A3;
    double ***Dr_th1_A4, ***Dr_th2_A4, ***Dr_ph_A4;
    double ***Dth_ph_A1, ***Dth_r_A1;
    double ***Dth_ph_A2, ***Dth_r_A2;
    double ***Dth_ph_A3, ***Dth_r_A3;
    double ***Dth_ph_A4, ***Dth_r_A4;
    double ***Dph_r_A1, ***Dph_th_A1;
    double ***Dph_r_A2, ***Dph_th_A2;
    double ***Dph_r_A3, ***Dph_th_A3;
    double ***Dph_r_A4, ***Dph_th_A4;

    double ***Hr_th1_A1, ***Hr_th2_A1, ***Hr_ph_A1;
    double ***Hr_th1_A2, ***Hr_th2_A2, ***Hr_ph_A2;
    double ***Hr_th1_A3, ***Hr_th2_A3, ***Hr_ph_A3;
    double ***Hr_th1_A4, ***Hr_th2_A4, ***Hr_ph_A4;
    double ***Hth_ph_A1, ***Hth_r_A1;
    double ***Hth_ph_A2, ***Hth_r_A2;
    double ***Hth_ph_A3, ***Hth_r_A3;
    double ***Hth_ph_A4, ***Hth_r_A4;
    double ***Hph_r_A1, ***Hph_th_A1;
    double ***Hph_r_A2, ***Hph_th_A2;*/

    pml *idx_Dr = new pml[4];
    pml *idx_Dth = new pml[4];
    pml *idx_Dph = new pml[4];
    pml *idx_Hr = new pml[4];
    pml *idx_Hth = new pml[4];
    pml *idx_Hph = new pml[4];

    PML_idx_initialize( idx_Dr, idx_Dth, idx_Dph,
                        idx_Hr, idx_Hth, idx_Hph );

    sigma_calc( sigma_th, sigma_ph, sigma_th_h, sigma_ph_h );

    t = Dt*0.0;

    std::cout << "Nr : " << Nr << " Ntheta : " << Ntheta
                << " Nphi : " << Nphi << "\n";
    
    std::chrono::system_clock::time_point start
        = std::chrono::system_clock::now();
    
    for( int n = 0; n < time_step + 1; n++ ){

        std::cout << "\r" << n << " / " << time_step;
        
        t = n*Dt;

        // Add Jth //
        J = -((t - t0)/sigma_t/sigma_t/delta_r/(dist(i_s + 0.5)*delta_theta)/(dist(i_s + 0.5)*delta_phi))
            *std::exp(-std::pow(t - t0, 2.0)/2.0/std::pow(sigma_t, 2.0));
        
        nEth[i_s][j_s][k_s] = nEth[i_s][j_s][k_s] + J;

        D_update( 
            nDr, nDth, nDph, oDr, oDth, oDph,
            Hr, Hth, Hph );
        
        D_update_pml(
            nDr, nDth, nDph, Hr, Hth, Hph,
            Dr_th1, Dr_th2, Dr_ph,
            Dth_ph, Dth_r,
            Dph_r, Dph_th,
            sigma_th, sigma_ph, idx_Dr, idx_Dth, idx_Dph );

        E_update(
            nEr, nEth, nEph, oEr, oEth, oEph,
            nDr, nDth, nDph, oDr, oDth, oDph );

        H_update(
            nEr, nEth, nEph, Hr, Hth, Hph );

        H_update_pml(
            nEr, nEth, nEph, Hr, Hth, Hph,
            Hr_th1, Hr_th2, Hr_ph, Hth_ph, Hth_r, Hph_r, Hph_th,
            sigma_th_h, sigma_ph_h, idx_Hr, idx_Hth, idx_Hph );

        std::string filename = "./result/ez_" << std::to_string(n) << ".dat";
        std::ofstream ofs(filename.c_str());
        for( int i = 0; i < Nr; i++ ){
            for( int k = 1; k < Nphi; k++ ){
                ofs << i << " " << k << nEth[i][j_s][k] << "\n";
            }
            ofs << "\n";
        }
        ofs.close();

    }

    std::cout << nEr[i_s][j_s][k_s] << "\n";

}
