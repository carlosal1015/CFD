#include "1D_BTCS.hh"

void Riemann_period()
{
  double x_l = 0;
  double x_r = 1.0;
  int nx = 200;
  double dx = (x_r - x_l) / nx;
  std::vector<double> x(nx, 0); // temperory array by RK3 scheme

  double t = 0.25;
  double dt = 0.0001;
  int nt = std::ceil(t / dt);
  dt = t / nt;
  std::vector<std::vector<double>> u(
      nt + 1, std::vector<double>(
                  nx, 0));       // one timestep = one row, FDM has nx+1 columns
  std::vector<double> ut(nx, 0); // temperory array by RK3 scheme
  std::vector<double> r(nx, 0);  // general spatial FD

  for (int i = 0; i < nx; i++) {
    x[i] = (i + 0.5) * dx;
    u[0][i] = std::sin(2.0 * Pi * x[i]);
  }

  for (int j = 1; j < nt + 1; j++) // one time step
  {
    rhs_RM(nx, dx, u[j - 1], r);
    for (int i = 1; i < nx - 1; i++) {
      ut[i] = u[j - 1][i] + dt * r[i];
    }
    rhs_RM(nx, dx, ut, r);
    for (int i = 1; i < nx - 1; i++) {
      ut[i] = 0.75 * u[j - 1][i] + 0.25 * ut[i] + dt / 4.0 * r[i];
    }
    rhs_RM(nx, dx, ut, r);
    for (int i = 1; i < nx - 1; i++) {
      u[j][i] = u[j - 1][i] / 3.0 + 2.0 / 3.0 * ut[i] + dt * 2.0 / 3.0 * r[i];
    }
  }

  std::ofstream outfile("Riemann_period.dat");
  if (outfile.is_open()) {
    for (int i = 0; i < nx; i++) {
      outfile << u[nt][i] << " ";
    }
    outfile << std::endl;
  }
  else {
    std::cerr << "Error: unable to open file for writing" << std::endl;
  }
  return;
}

void rhs_RM(int nx, double dx, std::vector<double> u, std::vector<double> &r)
{
  std::vector<double> uL(nx + 1, 0);
  std::vector<double> uR(nx + 1, 0);
  // left and right side fluxes at the interface
  std::vector<double> fL(nx + 1, 0);
  std::vector<double> fR(nx + 1, 0);
  std::vector<double> f(nx + 1, 0);

  uL = wenoL_FS(nx, u);
  uR = wenoR_FS(nx, u);

  for (int i = 0; i < nx + 1; i++) {
    fL[i] = 0.5 * uL[i] * uL[i];
    fR[i] = 0.5 * uR[i] * uR[i];
  }

  Riemann(nx, u, uL, uR, f, fL, fR);

  for (int i = 0; i < nx; i++) {
    r[i] = -(f[i + 1] - f[i]) / dx;
  }
  return;
}

void Riemann(int nx, std::vector<double> u, std::vector<double> uL,
             std::vector<double> uR, std::vector<double> &f,
             std::vector<double> fL, std::vector<double> fR)
{
  std::vector<double> c(nx + 1, 0);
  for (int i = 1; i < nx; i++) {
    c[i] = std::max(std::abs(u[i - 1]), std::abs(u[i]));
  }
  c[0] = std::max(std::abs(u[0]), std::abs(u[nx - 1]));
  c[nx] = std::max(std::abs(u[0]), std::abs(u[nx - 1]));

  // Interface fluxes(Rusanov)
  for (int i = 0; i < nx + 1; i++) {
    f[i] = 0.5 * (fR[i] + fL[i]) - 0.5 * c[i] * (uR[i] - uL[i]);
  }
}