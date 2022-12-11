#include "cross-correlation.h"

#include <cassert>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <vector>

// Compute the unsigned log2 of an unsigned value.
static inline size_t uilog2(size_t value) {
  size_t log2 = (value & (value - 1)) != 0;
  while (value >>= 1) ++log2;
  return log2;
}

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

static void FastFourierTransformImpl(const complex_vec_t &x,
                                     complex_vec_t *output,
                                     const bool inverse) {
  const size_t num_samples = x.size();
  assert(output->size() == num_samples);
  assert((num_samples & (num_samples - 1)) == 0 && "number of samples not a power of 2!");

  // Let's compute log2 of the number of samples.
  const unsigned log2_num_samples = log2(num_samples);
  std::fill(output->begin(), output->end(), 0);

  // We have to do a bit-reverse copy of the input in the output.
  // What does it mean? It means we have to take the index value starting from zero,
  // reverse its bits, and use that number as output index!
  // For instance, 1010 -> 0101.
  for (size_t i = 0; i < x.size(); ++i) {
    const size_t ri = reverse_bits<size_t>(i, log2_num_samples);
    (*output)[ri] = x[i];
  }

  using namespace std::complex_literals;
  const real_t sign = inverse ? 1.0 : -1.0;

  // We can work just on the output vector.
  // The fft allows us to iterate logn times.
  // https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm.
  // This is the iterative radix-2 FFT algorithm implemented using bit-reversal permutation.
  for (size_t s = 0; s < log2_num_samples; ++s) {
    // s-index represents size of the "fft leaf".
    // For instance, for 8 samples, the first iteration will take care of computing the sub-ffts
    // of size 2 and will generate then 8 intermediate values that in the next iteration will be used
    // to compute the sub-ffts of size 4! A good way of understand this is drawing the butterfly diagram.
    // https://www.youtube.com/watch?v=1mVbZLHLaf0&t=1810s.
    // That is why we can use directly the output vector to compute the final fft!
    const size_t m = 1 << (s + 1);

    // Weight to be multiplied for every fft.
    std::complex<real_t> w_m = std::polar<real_t>(1, sign * M_PI * 2 / m);

    // First iterate through every sub-fft. Index "k" is the starting
    // index of the output vector for the sub-fft m. For instance the
    // first iteration will be 0, 2, 4, the second will be 0, 4, 8, etc..
    for (size_t k = 0; k < num_samples; k += m) {
      std::complex<real_t> w = 1;
      for (size_t j = 0; j < (m / 2); ++j) {
        // Butterfly diagram "diagonal" edges, they are multiplied by w.
        const auto t = w * (*output)[k + j + m / 2];

        // Butterfly diagram "direct" edges. Not multiplied by anything.
        const auto u = (*output)[k + j];

        // Update the values for this couple of sub-ffts for
        // both the sub-coefficients.
        (*output)[k + j] = u + t;
        (*output)[k + j + m / 2] = u - t;

        // Prepare w for the next sub-coefficients!
        w = w * w_m;
      }
    }
  }
  if (inverse) {
    real_t float_size = static_cast<real_t>(output->size());
    for (size_t i = 0; i < output->size(); ++i) {
      (*output)[i] /= float_size;
    }
  }
}

void InvFFT(const complex_vec_t &d, complex_vec_t *output) {
  return FastFourierTransformImpl(d, output, true);
}

void FFT(const complex_vec_t &d, complex_vec_t *output) {
  return FastFourierTransformImpl(d, output, false);
}

// Our kernel should match our input plus some padding to have a linear
// convolution. We assume the size of x is a power of 2. To understand why we
// shift (to have padding in the center) see:
// https://dsp.stackexchange.com/questions/82273/why-to-pad-zeros-at-the-middle-of-sequence-instead-at-the-end-of-the-sequence
template <typename Iter>
static std::vector<std::complex<real_t>> Pad(Iter begin, Iter end, const bool shift) {
  const size_t size = end - begin;
  std::vector<std::complex<real_t>> out(2 * size, 0);
  const size_t offset = shift ? size + size / 2 : size / 2;
  for (unsigned i = 0; i < size; ++i) {
    out[(offset + i) % (2 * size)] = *begin++;
  }
  return out;
}

static bool PrintFirst(const char *msg) {
  fprintf(stderr, "%s sizeof(real_t)=%d\n", msg, (int)sizeof(real_t));
  return true;
}

void PreprocessCorrelate(const std::span<real_t> &data,
                         correlate_preprocessed_t *result) {
  if (result->first.empty()) {
    result->first.resize(2*data.size());
    result->second.resize(2*data.size());
  }
  // Pad and reverse. We are padding because we might
  // want to correlate with a shorter filter b.
  // We also reverse the filter because we are performing a convolution.

  const auto padded_d = Pad(data.begin(), data.end(), false);
  FFT(padded_d, &result->first);

  const auto padded_revesed_d = Pad(data.rbegin(), data.rend(), true);
  FFT(padded_revesed_d, &result->second);
}

const std::vector<real_t> cross_correlate(const correlate_preprocessed_t &a,
                                          const correlate_preprocessed_t &b) {
  static bool init = PrintFirst("FFT cross correlation using preprocess");
  const complex_vec_t &padded_fft_a = a.first;
  const complex_vec_t &padded_fft_b = b.second;

  complex_vec_t multiplaction_result(padded_fft_a.size());
  for (unsigned i = 0; i < padded_fft_a.size(); ++i) {
    multiplaction_result[i] = padded_fft_a[i] * padded_fft_b[i];
  }
  complex_vec_t c(padded_fft_a.size());
  InvFFT(multiplaction_result, &c);
  const size_t num_samples = padded_fft_a.size()/2;

  // TODO: do the following two things in one step.
  std::vector<real_t> convolved(num_samples);
  for (size_t i = 0; i < convolved.size(); ++i) {
    convolved[i] = c[i + num_samples / 2 - 1].real();
  }

  std::vector<real_t> result(num_samples);
  const int offset = convolved.size()/2;
  for (size_t i = 0; i < convolved.size()/2; ++i) {
    result[i] = convolved[offset + i];
  }
  return result;
}
