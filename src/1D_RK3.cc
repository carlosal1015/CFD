#include "1D_BTCS.hh"

void RK3_1D()
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
  std::vector<std::vector<double>> u_onethird(
      nt + 1, std::vector<double>(nx + 1, 0)); // one timestep = one row
  std::vector<std::vector<double>> u_twothird(
      nt + 1, std::vector<double>(nx + 1, 0)); // one timestep = one row

  // initial condition
  for (int i = 0; i < nx + 1; i++) {
    u[0][i] = -std::sin(Pi * (dx * i - 1));
  }
  for (int j = 1; j < nt + 1; j++) // one time step
  {
    u_onethird[j][0] = 0.0; // BC
    for (int i = 1; i < nx; i++) {
      u_onethird[j][i] =
          u[j - 1][i] +
          CFL / 3 * (u[j - 1][i + 1] - 2 * u[j - 1][i] + u[j - 1][i - 1]);
    }
    u_onethird[j][nx] = 0.0; // BC
    u_twothird[j][0] = 0.0;
    for (int i = 1; i < nx; i++) {
      u_twothird[j][i] =
          u[j - 1][i] + CFL * 2 / 3 *
                            (u_onethird[j][i + 1] - 2 * u_onethird[j][i] +
                             u_onethird[j][i - 1]);
    }
    u_twothird[j][nx] = 0.0;
    u[j][0] = 0.0;
    for (int i = 1; i < nx; i++) {
      u[j][i] = u[j - 1][i] +
                CFL / 4 *
                    ((u[j - 1][i + 1] - 2 * u[j - 1][i] + u[j - 1][i - 1]) +
                     3 * (u_twothird[j][i + 1] - 2 * u_twothird[j][i] +
                          u_twothird[j][i - 1]));
    }
    u[j][nx] = 0.0;
  }

  std::ofstream outfile("1d_RK3_u.dat");
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