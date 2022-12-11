
#ifndef CROSS_CORRELATATION_H
#define CROSS_CORRELATATION_H

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

typedef float real_t;
typedef std::vector<std::complex<real_t>> complex_vec_t;
typedef std::pair<complex_vec_t, complex_vec_t> correlate_preprocessed_t;

void FFT(const complex_vec_t &d, complex_vec_t *output);
void InvFFT(const complex_vec_t &d, complex_vec_t *output);

void PreprocessCorrelate(const std::span<real_t> &data,
                         correlate_preprocessed_t *result);

const std::vector<real_t> cross_correlate(const correlate_preprocessed_t &a,
                                          const correlate_preprocessed_t &b);

#endif  // CROSS_CORRELATATION_H
