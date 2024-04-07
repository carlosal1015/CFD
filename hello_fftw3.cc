// https://stackoverflow.com/q/78034800
#include <fftw3.h>
#include <iostream>

int main()
{
  // Size of input data
  const int N = 8;

  // Allocate input and output arrays
  double *in = (double *)fftw_malloc(sizeof(double) * N);
  fftw_complex *out =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));

  // Create a plan for forward FFT
  fftw_plan plan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);

  // Initialize input data (example)
  for (int i = 0; i < N; ++i) {
    in[i] = i; // Example input data
  }

  // Execute FFT
  fftw_execute(plan);

  // Output FFT result
  for (int i = 0; i < N / 2 + 1; ++i) {
    std::cout << "FFT[" << i << "] = " << out[i][0] << " + " << out[i][1] << "i"
              << std::endl;
  }

  // Clean up
  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);
}