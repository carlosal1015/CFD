#include "1D_BTCS.hh"

void RK2_1D()
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
  double CFL = alpha * dt / std::pow(dx, 2);
  std::vector<std::vector<double>> u(
      nt + 1, std::vector<double>(nx + 1, 0)); // one timestep = one row
  std::vector<std::vector<double>> u_half(
      nt + 1, std::vector<double>(nx + 1, 0)); // one timestep = one row

  // initial condition
  for (int i = 0; i < nx + 1; i++) {
    u[0][i] = -std::sin(Pi * (dx * i - 1));
  }
  for (int j = 1; j < nt + 1; j++) // one time step
  {
    u_half[j][0] = 0.0; // BC
    for (int i = 1; i < nx; i++) {
      u_half[j][i] =
          u[j - 1][i] +
          CFL / 2 * (u[j - 1][i + 1] - 2 * u[j - 1][i] + u[j - 1][i - 1]);
    }
    u_half[j][nx] = 0.0; // BC
    u[j][0] = 0.0;
    for (int i = 1; i < nx; i++) {
      u[j][i] = u[j - 1][i] +
                CFL * (u_half[j][i + 1] - 2 * u_half[j][i] + u_half[j][i - 1]);
    }
    u[j][nx] = 0.0;
  }

  std::ofstream outfile("1d_RK2_u.dat");
  if (outfile.is_open()) {
    /*std::vector<std::vector<double>>::iterator ptr;
    for (ptr = u.begin(); ptr < u.end(); ptr++)
    {
            std::ostream_iterator<double> output_iterator(outfile, "\n");
            std::copy((*ptr).begin(), (*ptr).end(), output_iterator);
    }
    outfile.close();*/
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