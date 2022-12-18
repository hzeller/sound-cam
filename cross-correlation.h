
#ifndef CROSS_CORRELATATION_H
#define CROSS_CORRELATATION_H

#include <complex>
#include <cstddef>
#include <span>
#include <vector>
#include <span>

typedef float real_t;
typedef std::vector<std::complex<real_t>> complex_vec_t;
typedef std::span<std::complex<real_t>> complex_span_t;

void FFT(const complex_span_t in, complex_vec_t *out);
void InvFFT(const complex_span_t in, complex_vec_t *out);

// Print a bunch of arrays for gnuplottability. X axis is index, multiple
// columns with the real() part.
void PrintArray(FILE *out, const std::initializer_list<complex_span_t> &what);

#endif  // CROSS_CORRELATATION_H
