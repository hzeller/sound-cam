#include "cross-correlation.h"

#include <cassert>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <vector>


fftw_plan InvFFT(const complex_span_t in, complex_vec_t *out) {
  assert(in.size() == out->size());
  fftw_complex *in_data = (fftw_complex*) in.data();
  fftw_complex *out_data = (fftw_complex*) out->data();
  return fftw_plan_dft_1d(in.size(), in_data, out_data,
                          FFTW_BACKWARD, FFTW_MEASURE);
}

fftw_plan FFT(const complex_span_t in, complex_vec_t *out) {
  assert(in.size() == out->size());
  fftw_complex *in_data = (fftw_complex*) in.data();
  fftw_complex *out_data = (fftw_complex*) out->data();
  return fftw_plan_dft_1d(in.size(), in_data, out_data,
                          FFTW_FORWARD, FFTW_MEASURE);
}

void PrintArray(FILE *out, const std::initializer_list<complex_span_t> &what) {
  bool index_used = true;
  for (size_t i = 0; index_used; ++i) {
    fprintf(out, "%d ", (int)i);
    index_used = false;
    for (const auto &s : what) {
      if (i >= s.size()) {
        fprintf(out, " -");
      }
      else {
        index_used = true;
        fprintf(out, " %.4f", real_part(s[i]));
      }
    }
    fprintf(out, "\n");
  }
}
