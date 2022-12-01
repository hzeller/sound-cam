
#ifndef CROSS_CORRELATATION_H
#define CROSS_CORRELATATION_H

#include <vector>
#include <cstddef>

typedef float real_t;

// Cross-correlate the first "elements" of a and b, return cross correlation
// with "output_count" elements (thus it is possible to limit the effort if
// only a few output elements are needed).
// if "output_count" < 0, then it is assumed to be "elements".
// Vector "a" needs to have at least elements + output_count elements.
std::vector<real_t> cross_correlate(const std::vector<real_t> &a,
                                    const std::vector<real_t> &b,
                                    size_t elements,
                                    int output_count = -1);

#endif  // CROSS_CORRELATATION_H
