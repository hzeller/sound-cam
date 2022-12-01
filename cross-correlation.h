
#ifndef CROSS_CORRELATATION_H
#define CROSS_CORRELATATION_H

#include <vector>
#include <cstddef>
#include <complex>

typedef double real_t;
typedef std::vector<std::complex<real_t>> complex_vec_t;

complex_vec_t InverseFastFourierTransform(const complex_vec_t &x);
complex_vec_t FastFourierTransform(const complex_vec_t &x);

typedef std::pair<complex_vec_t, complex_vec_t> correlate_preprocessed_t;
correlate_preprocessed_t PreprocessCorrelate(const std::vector<real_t> &data);

const std::vector<real_t> cross_correlate(const correlate_preprocessed_t &a,
                                          const correlate_preprocessed_t &b);

#endif  // CROSS_CORRELATATION_H
