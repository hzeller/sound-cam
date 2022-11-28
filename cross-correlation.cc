#include <cassert>
#include <cstddef>
#include <vector>

std::vector<float> cross_correlate(const std::vector<float> &a,
                                   const std::vector<float> &b,
                                   size_t elements,
                                   int output_count) {
  if (output_count < 0) output_count = elements;
  
  assert(output_count <= (int)elements);
  assert(a.size() >= elements + output_count);
  assert(b.size() >= elements);

  std::vector<float> result(output_count);
  
  for (int i = 0; i < output_count; ++i) {
    for (size_t j = 0; j < elements; ++j) {
      result[i] += a[j + i] * b[j];
    }
  }
  return result;
}
