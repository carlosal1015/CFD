cmake_minimum_required(VERSION 3.20.0)
project("CFD"
  VERSION 1.0
  LANGUAGES CXX)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are not allowed")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(FTCS FTCS.cc
  src/1D_FTCS.cc)
target_include_directories(FTCS PUBLIC "include")

add_executable(BTCS BTCS.cc
  src/1D_BTCS.cc
  src/algebraicEqs.cc)
target_include_directories(BTCS PUBLIC "include")

add_executable(CNTCS CNTCS.cc
  src/1D_CrankNicolson.cc
  src/algebraicEqs.cc)
target_include_directories(CNTCS PUBLIC "include")

add_executable(RK2 RK2.cc
  src/1D_RK2.cc)
target_include_directories(RK2 PUBLIC "include")

add_executable(RK3 RK3.cc
  src/1D_RK3.cc)
target_include_directories(RK3 PUBLIC "include")

add_executable(RK4 RK4.cc
  src/1D_RK4.cc)
target_include_directories(RK4 PUBLIC "include")

add_executable(WENO5 WENO5.cc
  src/WENO5_Dirichlet.cc)
target_include_directories(WENO5 PUBLIC "include")

add_executable(CRWENO5 CRWENO5.cc
  src/Burgers_CRWENO_Dirichlet.cc)
target_include_directories(CRWENO5 PUBLIC "include")

add_executable(FluxSplitting FluxSplitting.cc
  src/Burgers_FluxSplitting_period.cc
  src/Burgers_WENO5_Dirichlet.cc)
target_include_directories(FluxSplitting PUBLIC "include")

add_executable(Riemann Riemann.cc
  src/Burgers_Riemann_period.cc
  src/Burgers_FluxSplitting_period.cc
  src/Burgers_WENO5_Dirichlet.cc)
target_include_directories(Riemann PUBLIC "include")

add_executable(Euler_roe Euler_roe.cc
  src/Euler_roe.cc
  src/Burgers_WENO5_Dirichlet.cc)
target_include_directories(Euler_roe PUBLIC "include")

add_executable(Euler_hllc Euler_hllc.cc
  src/Euler_HLLC.cc
  src/Euler_roe.cc
  src/Burgers_WENO5_Dirichlet.cc)
target_include_directories(Euler_hllc PUBLIC "include")

add_executable(Euler_rusanov Euler_rusanov.cc
  src/Euler_rusanov.cc
  src/Euler_roe.cc
  src/Burgers_WENO5_Dirichlet.cc)
target_include_directories(Euler_rusanov PUBLIC "include")

add_executable(hello_fftw3 hello_fftw3.cc)
target_link_libraries(hello_fftw3 PUBLIC fftw3)

add_executable(Poisson_FFT Poisson_FFT.cc
  src/Poisson_FFT.cc)
target_link_libraries(Poisson_FFT PUBLIC fftw3)
target_include_directories(Poisson_FFT PUBLIC "include")

add_executable(Poisson_CG Poisson_CG.cc
  src/Poisson_CG.cc)
target_link_libraries(Poisson_CG PUBLIC fftw3)
target_include_directories(Poisson_CG PUBLIC "include")

add_executable(Poisson_GaussSeidel Poisson_GaussSeidel.cc
  src/Poisson_GaussSeidel.cc)
target_link_libraries(Poisson_GaussSeidel PUBLIC fftw3)
target_include_directories(Poisson_GaussSeidel PUBLIC "include")

add_executable(Poisson_Multigrid Poisson_Multigrid.cc
  src/Poisson_Multigrid.cc)
target_link_libraries(Poisson_Multigrid PUBLIC fftw3)
target_include_directories(Poisson_Multigrid PUBLIC "include")