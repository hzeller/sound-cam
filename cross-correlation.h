
#ifndef CROSS_CORRELATATION_H
#define CROSS_CORRELATATION_H

#include <vector>
#include <cstddef>
#include <complex>

typedef double real_t;
typedef std::vector<std::complex<real_t>> complex_vec_t;

// Cross-correlate the first "elements" of a and b, return cross correlation
// with "output_count" elements (thus it is possible to limit the effort if
// only a few output elements are needed).
// if "output_count" < 0, then it is assumed to be "elements".
// Vector "a" needs to have at least elements + output_count elements.
std::vector<real_t> cross_correlate(const std::vector<real_t> &a,
                                    const std::vector<real_t> &b,
                                    size_t elements,
                                    int output_count = -1);

complex_vec_t InverseFastFourierTransform(const complex_vec_t &x);
complex_vec_t FastFourierTransform(const complex_vec_t &x);

typedef std::pair<complex_vec_t, complex_vec_t> correlate_preprocessed_t;
correlate_preprocessed_t PreprocessCorrelate(const std::vector<real_t> &data);

const std::vector<real_t> cross_correlate(const correlate_preprocessed_t &a,
                                          const correlate_preprocessed_t &b,
                                          int output_count);

#endif  // CROSS_CORRELATATION_H
