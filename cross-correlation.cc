#include <cassert>
#include <cstddef>
#include <vector>

#include "cross-correlation.h"

#include <alglib/fasttransforms.h>

static bool PrintFirst(const char *msg) {
  fprintf(stderr, "%s sizeof(real_t)=%d\n", msg, (int)sizeof(real_t));
  return true;
}

// Very simplistic O(N²) implementation.
#if CROSS_MANUAL_IMPL
std::vector<real_t> cross_correlate(const std::vector<real_t> &a,
                                    const std::vector<real_t> &b,
                                    size_t elements,
                                    int output_count) {
  static bool init = PrintFirst("Manual cross correlation");
  if (output_count < 0) output_count = elements;

  assert(output_count <= (int)elements);
  assert(a.size() >= elements + output_count);
  assert(b.size() >= elements);

  std::vector<real_t> result(output_count);

  for (int i = 0; i < output_count; ++i) {
    for (size_t j = 0; j < elements; ++j) {
      result[i] += a[j + i] * b[j];
    }
  }
  return result;
}
#else
std::vector<real_t> cross_correlate(const std::vector<real_t> &a,
                                    const std::vector<real_t> &b,
                                    size_t elements,
                                    int output_count) {
  static bool init = PrintFirst("alglib cross correlation");
  if (output_count < 0) output_count = elements;

  alglib::real_1d_array out;
  alglib::real_1d_array signal, pattern;
  signal.attach_to_ptr(a.size(), (real_t*)a.data());
  pattern.attach_to_ptr(b.size(), (real_t*)b.data());

  alglib::corrr1d(signal, signal.length(), pattern, elements, out);

  std::vector<real_t> result(output_count);
  std::copy(out.getcontent(), out.getcontent() + output_count,
            result.begin());
  return result;
}
#endif
