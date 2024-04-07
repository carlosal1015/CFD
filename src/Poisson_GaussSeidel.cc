#include "1D_BTCS.hh"
#include <complex>
#include <fftw3.h>

void Poisson_GaussSeidel()
{
  double x_l = 0.0;
  double x_r = 1.0;
  int nx = 512;
  double dx = (x_r - x_l) / nx;
  std::vector<double> x(nx + 1, 0);
  for (int i = 0; i < nx + 1; i++) {
    x[i] = i * dx + x_l;
  }

  double y_b = 0.0;
  double y_t = 1.0;
  int ny = 512;
  double dy = (y_t - y_b) / ny;
  std::vector<double> y(ny + 1, 0);
  for (int i = 0; i < ny + 1; i++) {
    y[i] = i * dy + y_b;
  }

  double tolerance = 1.0e-4;
  int max_iter = 10000;

  std::vector<std::vector<double>> ue(ny + 1, std::vector<double>(nx + 1, 0.0));
  std::vector<std::vector<double>> f(ny + 1, std::vector<double>(nx + 1, 0.0));
  std::vector<std::vector<double>> un(ny + 1, std::vector<double>(nx + 1, 0.0));

  // analytic solution and initial condition
  for (int i = 0; i < ny + 1; i++) {
    for (int j = 0; j < nx + 1; j++) {
      ue[i][j] = (x[j] * x[j] - 1.0) * (y[i] * y[i] - 1.0);
      f[i][j] = -2.0 * (2.0 - x[j] * x[j] - y[i] * y[i]);
    }
  }
  for (int i = 0; i < ny + 1; i++) {
    un[i][0] = ue[i][0];
    un[i][ny] = ue[i][ny];
  }
  for (int i = 0; i < nx + 1; i++) {
    un[0][i] = ue[0][i];
    un[nx][i] = ue[nx][i];
  }

  std::vector<std::vector<double>> r(ny + 1, std::vector<double>(nx + 1, 0.0));
  double init_rms = 0.0;
  double rms = 0.0;

  for (int i = 1; i < ny; i++) {
    for (int j = 1; j < nx; j++) {
      double d2udx2 = (un[i + 1][j] - 2 * un[i][j] + un[i - 1][j]) / dx / dx;
      double d2udy2 = (un[i][j + 1] - 2 * un[i][j] + un[i][j - 1]) / dy / dy;
      r[i][j] = f[i][j] - d2udx2 - d2udy2;
    }
  }
  // compute residual
  for (int i = 1; i < ny; i++) {
    for (int j = 1; j < nx; j++) {
      init_rms += r[i][j] * r[i][j];
    }
  }
  init_rms = std::sqrt(init_rms / (nx - 1) / (ny - 1));
  rms = init_rms;

  int iter_count = 0;
  double den = -2.0 / dx / dx - 2.0 / dy / dy;
  double exp_rms = tolerance * init_rms;

  for (iter_count = 0; iter_count < max_iter && rms > exp_rms; iter_count++) {
    // correct solution
    for (int i = 1; i < ny; i++) {
      for (int j = 1; j < nx; j++) {
        double d2udx2 = (un[i + 1][j] - 2 * un[i][j] + un[i - 1][j]) / dx / dx;
        double d2udy2 = (un[i][j + 1] - 2 * un[i][j] + un[i][j - 1]) / dy / dy;
        r[i][j] = f[i][j] - d2udx2 - d2udy2;
        un[i][j] += r[i][j] / den;
      }
    }
    // compute new residual
    rms = 0.0;
    for (int i = 1; i < ny; i++) {
      for (int j = 1; j < nx; j++) {
        double d2udx2 = (un[i + 1][j] - 2 * un[i][j] + un[i - 1][j]) / dx / dx;
        double d2udy2 = (un[i][j + 1] - 2 * un[i][j] + un[i][j - 1]) / dy / dy;
        r[i][j] = f[i][j] - d2udx2 - d2udy2;
        rms += r[i][j] * r[i][j];
      }
    }
    rms = std::sqrt(rms / (nx - 1) / (ny - 1));
    std::cout << "iteration times " << iter_count << std::endl;
    std::cout << "residual " << rms << std::endl;
  }
  std::cout << "iteration times until convergence:" << iter_count << std::endl;

  // write
  std::ofstream outfile("Poisson_GaussSeidel.dat");
  if (outfile.is_open()) {
    for (int i = 0; i < ny + 1; i++) {
      for (int j = 0; j < nx + 1; j++) {
        outfile << un[i][j] << " ";
      }
      outfile << std::endl;
    }
    outfile << std::endl;
  }
  else {
    std::cerr << "Error: unable to open file for writing" << std::endl;
  }
  return;
}
