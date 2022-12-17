
#ifndef CROSS_CORRELATATION_H
#define CROSS_CORRELATATION_H

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

typedef float real_t;
typedef std::vector<std::complex<real_t>> complex_vec_t;

void FFT(const complex_vec_t &in, complex_vec_t *out);
void InvFFT(const complex_vec_t &in, complex_vec_t *out);

#endif  // CROSS_CORRELATATION_H
