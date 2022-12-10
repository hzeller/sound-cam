
#ifndef CROSS_CORRELATATION_H
#define CROSS_CORRELATATION_H

#include <vector>
#include <cstddef>
#include <complex>

typedef float real_t;
typedef std::vector<std::complex<real_t>> complex_vec_t;
typedef std::pair<complex_vec_t, complex_vec_t> correlate_preprocessed_t;

void PreprocessCorrelate(const std::vector<real_t> &data,
                         correlate_preprocessed_t *result);

const std::vector<real_t> cross_correlate(const correlate_preprocessed_t &a,
                                          const correlate_preprocessed_t &b);

#endif  // CROSS_CORRELATATION_H
