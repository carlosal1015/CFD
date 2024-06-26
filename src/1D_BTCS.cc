#include "1D_BTCS.hh"

void BTCS_1D()
{
  double x_l = -1.0;
  double x_r = 1.0;
  double dx = 0.025;
  int nx = std::ceil((x_r - x_l) / dx);
  dx = (x_r - x_l) / nx;

  double t = 1.0;
  double dt = 0.0025;
  int nt = std::ceil(t / dt);
  dt = t / nt;

  double alpha = 1 / (Pi * Pi);
  double CFL = alpha * dt / pow(dx, 2);
  std::vector<std::vector<double>> u(
      nt + 1, std::vector<double>(nx + 1, 0)); // one timestep = one row
  std::vector<double> a(nx + 1, 0);
  std::vector<double> b(nx + 1, 0);
  std::vector<double> c(nx + 1, 0);
  std::vector<double> q(nx + 1, 0);

  // initial condition
  for (int i = 0; i < nx + 1; i++) {
    u[0][i] = -std::sin(Pi * (dx * i - 1));
  }
  a[0] = c[0] = q[0] = 0;
  b[0] = 1;
  for (int i = 1; i < nx; i++) {
    a[i] = c[i] = -CFL;
    b[i] = 1 + 2 * CFL;
  }
  a[nx] = c[nx] = q[nx] = 0;
  b[nx] = 1;
  for (int j = 1; j < nt + 1; j++) {
    for (int i = 1; i < nx; i++) {
      q[i] = u[j - 1][i];
    }
    thomasTridiagonal(a, b, c, u[j], q);
  }

  std::ofstream outfile("1d_BTCS_u.dat");
  if (outfile.is_open()) {
    for (int i = 0; i < nx + 1; i++) {
      outfile << u[nt][i] << " ";
    }
    outfile << std::endl;
  }
  else {
    std::cerr << "Error: unable to open file for writing" << std::endl;
  }
  return;
}