#include "1D_BTCS.hh"

void FluxSplitting_Burgers()
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
    u[j][0] = 0.0;      // BC
    u[j][nx - 1] = 0.0; // BC

    rhs_FS(nx, dx, u[j - 1], r);
    for (int i = 1; i < nx - 1; i++) {
      ut[i] = u[j - 1][i] + dt * r[i];
    }
    rhs_FS(nx, dx, ut, r);
    for (int i = 1; i < nx - 1; i++) {
      ut[i] = 0.75 * u[j - 1][i] + 0.25 * ut[i] + dt / 4.0 * r[i];
    }
    rhs_FS(nx, dx, ut, r);
    for (int i = 1; i < nx - 1; i++) {
      u[j][i] = u[j - 1][i] / 3.0 + 2.0 / 3.0 * ut[i] + dt * 2.0 / 3.0 * r[i];
    }
  }

  std::ofstream outfile("FluxSplitting_Burgers.dat");
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

void rhs_FS(int nx, double dx, std::vector<double> u, std::vector<double> &r)
{
  // flux computed at nodal pointsand positiveand negative splitting
  std::vector<double> f(nx, 0);
  std::vector<double> fP(nx, 0); // flux positive
  std::vector<double> fN(nx, 0);

  // wave speed at nodal points
  std::vector<double> alpha(nx, 0);

  // left and right side fluxes at the interface
  std::vector<double> fL(nx + 1, 0);
  std::vector<double> fR(nx + 1, 0);

  for (int i = 0; i < nx; i++) {
    f[i] = 0.5 * u[i] * u[i];
  }
  waveSpeed(nx, u, alpha); // df/du
  for (int i = 0; i < nx; i++) {
    fP[i] = 0.5 * (f[i] + alpha[i] * u[i]);
    fN[i] = 0.5 * (f[i] - alpha[i] * u[i]);
  }
  fL = wenoL_FS(nx, fP);
  fR = wenoR_FS(nx, fN);

  for (int i = 0; i < nx; i++) {
    r[i] = -(fL[i + 1] - fL[i]) / dx - (fR[i + 1] - fR[i]) / dx;
  }
  return;
}

void waveSpeed(int nx, std::vector<double> u, std::vector<double> &alpha)
{
  int i = 2;
  for (i = 2; i < nx - 2; i++) {
    alpha[i] =
        std::max(std::max(std::max(std::abs(u[i - 2]), std::abs(u[i - 1])),
                          std::max(std::abs(u[i]), std::abs(u[i + 1]))),
                 std::abs(u[i + 2]));
  }
  // periodic BC
  i = 0;
  alpha[i] =
      std::max(std::max(std::max(std::abs(u[nx - 2]), std::abs(u[nx - 1])),
                        std::max(std::abs(u[i]), std::abs(u[i + 1]))),
               std::abs(u[i + 2]));
  i = 1;
  alpha[i] =
      std::max(std::max(std::max(std::abs(u[nx - 1]), std::abs(u[i - 1])),
                        std::max(std::abs(u[i]), std::abs(u[i + 1]))),
               std::abs(u[i + 2]));
  i = nx - 2;
  alpha[i] = std::max(std::max(std::max(std::abs(u[i - 2]), std::abs(u[i - 1])),
                               std::max(std::abs(u[i]), std::abs(u[i + 1]))),
                      std::abs(u[0]));
  i = nx - 1;
  alpha[i] = std::max(std::max(std::max(std::abs(u[i - 2]), std::abs(u[i - 1])),
                               std::max(std::abs(u[i]), std::abs(u[0]))),
                      std::abs(u[1]));
}
std::vector<double> wenoL_FS(int nx, std::vector<double> u)
{
  std::vector<double> f(nx + 1, 0);

  int i = -1;
  double v1 = u[nx - 3];
  double v2 = u[nx - 2];
  double v3 = u[nx - 1];
  double v4 = u[i + 1];
  double v5 = u[i + 2];
  f[i + 1] = wcL(v1, v2, v3, v4, v5);

  i = 0;
  v1 = u[nx - 2];
  v2 = u[nx - 1];
  v3 = u[i];
  v4 = u[i + 1];
  v5 = u[i + 2];
  f[i + 1] = wcL(v1, v2, v3, v4, v5);

  i = 1;
  v1 = u[nx - 1];
  v2 = u[i - 1];
  v3 = u[i];
  v4 = u[i + 1];
  v5 = u[i + 2];
  f[i + 1] = wcL(v1, v2, v3, v4, v5);

  for (i = 2; i < nx - 2; i++) {
    v1 = u[i - 2];
    v2 = u[i - 1];
    v3 = u[i];
    v4 = u[i + 1];
    v5 = u[i + 2];
    f[i + 1] = wcL(v1, v2, v3, v4, v5);
  }
  i = nx - 2;
  v1 = u[i - 2];
  v2 = u[i - 1];
  v3 = u[i];
  v4 = u[i + 1];
  v5 = u[0];
  f[i + 1] = wcL(v1, v2, v3, v4, v5);

  i = nx - 1;
  v1 = u[i - 2];
  v2 = u[i - 1];
  v3 = u[i];
  v4 = u[0];
  v5 = u[1];
  f[i + 1] = wcL(v1, v2, v3, v4, v5);

  return f;
}

std::vector<double> wenoR_FS(int nx, std::vector<double> u)
{
  std::vector<double> f(nx + 1, 0);
  int i = 0;
  double v1 = u[nx - 2];
  double v2 = u[nx - 1];
  double v3 = u[i];
  double v4 = u[i + 1];
  double v5 = u[i + 2];
  f[i] = wcR(v1, v2, v3, v4, v5);

  i = 1;
  v1 = u[nx - 1];
  v2 = u[i - 1];
  v3 = u[i];
  v4 = u[i + 1];
  v5 = u[i + 2];
  f[i] = wcR(v1, v2, v3, v4, v5);

  for (i = 2; i < nx - 2; i++) {
    v1 = u[i - 2];
    v2 = u[i - 1];
    v3 = u[i];
    v4 = u[i + 1];
    v5 = u[i + 2];
    f[i] = wcR(v1, v2, v3, v4, v5);
  }

  i = nx - 2;
  v1 = u[i - 2];
  v2 = u[i - 1];
  v3 = u[i];
  v4 = u[i + 1];
  v5 = u[0];
  f[i] = wcR(v1, v2, v3, v4, v5);

  i = nx - 1;
  v1 = u[i - 2];
  v2 = u[i - 1];
  v3 = u[i];
  v4 = u[0];
  v5 = u[1];
  f[i] = wcR(v1, v2, v3, v4, v5);

  i = nx;
  v1 = u[i - 2];
  v2 = u[i - 1];
  v3 = u[0];
  v4 = u[1];
  v5 = u[2];
  f[i] = wcR(v1, v2, v3, v4, v5);

  return f;
}