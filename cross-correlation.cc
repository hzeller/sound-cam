#include "cross-correlation.h"

#include <cassert>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <vector>

// Naive approach to reverse the first <width> LSBs of an integer n.
template <class T>
T reverse_bits(T n, short width) {
  T output = T(0);
  for (int i = 0; i < width; ++i) {
    const short pos = width - i - 1;
    output |= ((n >> pos) & 1) << i;
  }
  return output;
}

static void FastFourierTransformImpl(complex_span_t in,
                                     complex_vec_t *out,
                                     const bool inverse) {
  const size_t num_samples = in.size();
  assert(out->size() == num_samples);
  assert((num_samples & (num_samples - 1)) == 0 &&
         "number of samples not a power of 2!");

  // Let's compute log2 of the number of samples.
  const unsigned log2_num_samples = log2(num_samples);

  // We have to do a bit-reverse copy of the input in the output.
  // What does it mean? It means we have to take the index value starting
  // from zero, reverse its bits, and use that number as output index!
  // For instance, 1010 -> 0101.
  for (size_t i = 0; i < num_samples; ++i) {
    const size_t ri = reverse_bits<size_t>(i, log2_num_samples);
    (*out)[ri] = in[i];
  }

  const real_t signed_tau = inverse ? 2.0 * M_PI : -2.0 * M_PI;

  // We can work just on the output vector.
  // The fft allows us to iterate logn times.
  // https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm.
  // This is the iterative radix-2 FFT algorithm implemented using
  // bit-reversal permutation.
  for (size_t s = 0; s < log2_num_samples; ++s) {
    // s-index represents size of the "fft leaf".
    // For instance, for 8 samples, the first iteration will take care of
    // computing the sub-ffts of size 2 and will generate then 8 intermediate
    // values that in the next iteration will be used to compute the
    // sub-ffts of size 4! A good way of understand this is drawing the
    // butterfly diagram.
    // https://www.youtube.com/watch?v=1mVbZLHLaf0&t=1810s
    // That is why we can use directly the output vector to compute
    // the final fft!
    const size_t m = 1 << (s + 1);

    // Weight to be multiplied for every fft.
    std::complex<real_t> w_m = std::polar<real_t>(1, signed_tau / m);

    // First iterate through every sub-fft. Index "k" is the starting
    // index of the output vector for the sub-fft m. For instance the
    // first iteration will be 0, 2, 4, the second will be 0, 4, 8, etc..
    for (size_t k = 0; k < num_samples; k += m) {
      std::complex<real_t> w = 1;
      for (size_t j = 0; j < (m / 2); ++j) {
        // Butterfly diagram "diagonal" edges, they are multiplied by w.
        const auto t = w * (*out)[k + j + m / 2];

        // Butterfly diagram "direct" edges. Not multiplied by anything.
        const auto u = (*out)[k + j];

        // Update the values for this couple of sub-ffts for
        // both the sub-coefficients.
        (*out)[k + j] = u + t;
        (*out)[k + j + m / 2] = u - t;

        // Prepare w for the next sub-coefficients!
        w = w * w_m;
      }
    }
  }
  if (inverse) {
    const real_t float_size = static_cast<real_t>(out->size());
    for (size_t i = 0; i < num_samples; ++i) {
      (*out)[i] /= float_size;
    }
  }
}

void InvFFT(const complex_span_t in, complex_vec_t *out) {
  return FastFourierTransformImpl(in, out, true);
}

void FFT(const complex_span_t in, complex_vec_t *out) {
  return FastFourierTransformImpl(in, out, false);
}
