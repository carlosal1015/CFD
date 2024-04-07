#pragma once
constexpr auto Pi = 3.1415926;

#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <ostream>
#include <vector>

// for heat transfer
void BTCS_1D();
void FTCS_1D();
void CNTCS_1D();
void RK2_1D();
void RK3_1D();
void RK4_1D();

// for solution of algebraic eqs.
void thomasTridiagonal(std::vector<double>, std::vector<double>,
                       std::vector<double>, std::vector<double> &,
                       std::vector<double>);

// Burger's equation
//  for WENO5
void WENO5_Dirichlet();
void rhs(int, double, std::vector<double>, std::vector<double> &);
void wenoL(int, std::vector<double>, std::vector<double> &);
void wenoR(int, std::vector<double>, std::vector<double> &);
double wcL(double, double, double, double, double);
double wcR(double, double, double, double, double);
// for CRWENO5
void CRWENO5_Dirichlet();
void rhs_CR(int, double, std::vector<double>,
            std::vector<double> &); // Calculate right hand term of the inviscid
                                    // Burgers equation
void wenoL_CR(int, std::vector<double>,
              std::vector<double> &); // CRWENO reconstruction ofr upwind
                                      // direction (positive and left to right)
void wenoR_CR(int, std::vector<double>, std::vector<double> &);
std::vector<double> wcL_CR(double, double, double, double,
                           double); // nonlinear weights for upwind direction
std::vector<double> wcR_CR(double, double, double, double, double);
void thomasTridiagonal_CRWENO(std::vector<double>, std::vector<double>,
                              std::vector<double>, std::vector<double> &,
                              std::vector<double>, int, int);
// for FVM-Flux splitting
void FluxSplitting_Burgers();
void rhs_FS(int, double, std::vector<double>, std::vector<double> &);
void waveSpeed(int, std::vector<double>, std::vector<double> &);
std::vector<double> wenoL_FS(int, std::vector<double>);
std::vector<double> wenoR_FS(int, std::vector<double>);
// for FVM-Riemann solver
void Riemann_period();
void rhs_RM(int, double, std::vector<double>, std::vector<double> &);
void Riemann(int, std::vector<double>, std::vector<double>, std::vector<double>,
             std::vector<double> &, std::vector<double>, std::vector<double>);

// 1D Euler Solver
void Euler_roe();
void rhs_roe(int, double, double, std::vector<std::vector<double>>,
             std::vector<std::vector<double>> &);
void flux_roe(int, double, std::vector<std::vector<double>>,
              std::vector<std::vector<double>> &);
void roe(int, double, std::vector<std::vector<double>>,
         std::vector<std::vector<double>>, std::vector<std::vector<double>> &,
         std::vector<std::vector<double>>, std::vector<std::vector<double>>);
std::vector<std::vector<double>> wenoL_roe(int,
                                           std::vector<std::vector<double>>);
std::vector<std::vector<double>> wenoR_roe(int,
                                           std::vector<std::vector<double>>);
// HLLC
void Euler_hllc();
void rhs_hllc(int, double, double, std::vector<std::vector<double>>,
              std::vector<std::vector<double>> &);
void hllc(int, double, std::vector<std::vector<double>>,
          std::vector<std::vector<double>>, std::vector<std::vector<double>> &,
          std::vector<std::vector<double>>, std::vector<std::vector<double>>);
// Rusanov
// HLLC
void Euler_rusanov();
void rhs_rusanov(int, double, double, std::vector<std::vector<double>>,
                 std::vector<std::vector<double>> &);
void rusanov(int, double, std::vector<std::vector<double>>,
             std::vector<std::vector<double>>,
             std::vector<std::vector<double>> &,
             std::vector<std::vector<double>>,
             std::vector<std::vector<double>>);

// 2D Poisson equation
void Poisson_FFT();
void Poisson_CG();