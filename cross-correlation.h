
#ifndef CROSS_CORRELATATION_H
#define CROSS_CORRELATATION_H

#include <vector>
#include <cstddef>

// Cross-correlate the first "elements" of a and b, return cross correlation
// with "output_count" elements (thus it is possible to limit the effort if
// only a few output elements are needed).
// if "output_count" < 0, then it is assumed to be "elements".
// Vector "a" needs to have at least elements + output_count elements.
std::vector<float> cross_correlate(const std::vector<float> &a,
                                   const std::vector<float> &b,
                                   size_t elements,
				   int output_count = -1);

#endif  // CROSS_CORRELATATION_H
