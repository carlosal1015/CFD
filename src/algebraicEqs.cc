#include "1D_BTCS.hh"

void thomasTridiagonal(std::vector<double> a, std::vector<double> b,
                       std::vector<double> c, std::vector<double> &u,
                       std::vector<double> q)
{
  int n_eq = a.size();
  std::vector<double> m(n_eq, 0);
  u[0] = q[0] / b[0];

  for (int i = 1; i < n_eq; i++) {
    m[i] = c[i - 1] / b[i - 1];
    b[i] = b[i] - m[i] * a[i];
    u[i] = (q[i] - u[i - 1] * a[i]) / b[i];
  }
  u[n_eq - 1] = q[n_eq - 1];
  for (int i = n_eq - 2; i >= 0; i--) {
    u[i] = u[i] - u[i + 1] * m[i + 1];
  }
  return;
}