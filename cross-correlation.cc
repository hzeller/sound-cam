#include "cross-correlation.h"

#include <cassert>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <vector>


fftwf_plan InvFFT(const complex_span_t in, complex_span_t *out) {
  assert(in.size() == out->size());
  // fftw_complex and std::complex<> have the same memory layout.
  fftwf_complex *in_data = (fftwf_complex*) in.data();
  fftwf_complex *out_data = (fftwf_complex*) out->data();
  return fftwf_plan_dft_1d(in.size(), in_data, out_data,
                          FFTW_BACKWARD, FFTW_MEASURE);
}

fftwf_plan InvFFTReal(const complex_span_t in, real_span_t *out) {
  assert(in.size() == out->size());
  // fftw_complex and std::complex<> have the same memory layout.
  fftwf_complex *in_data = (fftwf_complex*) in.data();
  real_t *out_data = (real_t*) out->data();
  return fftwf_plan_dft_c2r_1d(in.size(), in_data, out_data, FFTW_MEASURE);
}

fftwf_plan FFT(const complex_span_t in, complex_span_t *out) {
  assert(in.size() == out->size());
  fftwf_complex *in_data = (fftwf_complex*) in.data();
  fftwf_complex *out_data = (fftwf_complex*) out->data();
  return fftwf_plan_dft_1d(in.size(), in_data, out_data,
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
        fprintf(out, " %.4f", s[i].real());
      }
    }
    fprintf(out, "\n");
  }
}
