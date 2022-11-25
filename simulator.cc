/*
 * Simulation of sound cam
 */

#include <assert.h>
#include <chrono>
#include <functional>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <vector>

#include "colormap.h"
#include "terminal-canvas.h"

constexpr float tau = 2 * M_PI;
constexpr float rad2deg = 360 / tau;

constexpr float kSpeedOfSound = 340.0f; // m/s
constexpr int kNoiseBits = 2;
#if 1
// Regular
constexpr int kSampleBits = 12;
constexpr int kSampleRateHz = 48000; // Hz
#else
// high time resolution.
constexpr int kSampleBits = 12;
constexpr int kSampleRateHz = 192000; // Hz
#endif

constexpr float kMeasureTime = 1;              // second;
constexpr float kTestSourceFrequency = 1200.0; // baseline sound frequency.

constexpr float display_range = tau / 4; // Angle of view.
constexpr int resolution_pixels = 121;

#define cm *0.01f

static float square(float a) { return a * a; }
class Point {
public:
  static Point FromAngle(float yaw, float pitch) {
    return {cos(yaw) * cos(pitch), sin(yaw) * cos(pitch), sin(pitch)};
  }

  float distance_to(const Point &other) const {
    return sqrtf(square(x - other.x) + square(y - other.y) +
                 square(z - other.z));
  }

  float length() const { return sqrtf(square(x) + square(y) + square(z)); }

  float dotMul(const Point &other) const {
    return x * other.x + y * other.y + z * other.z;
  }

  void MakeUnitLen() {
    float len = length();
    x /= len;
    y /= len;
    z /= len;
  }

  float x, y, z;
};

Point operator+(const Point &a, const Point &b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Point operator-(const Point &a, const Point &b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Point operator*(const Point &p, float scalar) {
  return {p.x * scalar, p.y * scalar, p.z * scalar};
}

Point sound_source1 = Point({0.3, 0.3, 1}) * 1;
Point sound_source2 = Point({0, -0.4, 1}) * 1;

Point optical_camera_pos = {0, 0, 0};

#define arraysize(a) sizeof(a) / sizeof(a[0])

// Sound source. Returns samples from -1 to 1 at given time.
typedef std::function<float(float)> SoundSource;

class MicrophoneSamples {
public:
  MicrophoneSamples(const Point &microphone, int samples)
      : mic_pos_(microphone), values_(samples) {}

  const Point &pos() const { return mic_pos_; }

  void fill(SoundSource source, float distance_m) {
    const float sample_distance = 1.0f / kSampleRateHz;
    const float time_delay = distance_m / kSpeedOfSound;
    // We keep things symmetric around 0, so half above, half below; -1bit
    const int noisefloor_magnitude = (1 << kNoiseBits);
    // To keep things simple, we limit the range of the signal so that
    // it never overflows and give the noisefloor its own space.
    const int quant = (1 << (kSampleBits - kNoiseBits - 1));
    for (size_t i = 0; i < values_.size(); ++i) {
      const int value = quant * source(i * sample_distance - time_delay) +
                        random() % noisefloor_magnitude;
      values_[i] += (1.0 * value / quant);
    }
  }

  float value_at(float t) {
    // Simplest: find closest to sample point
    int sample = round(t * kSampleRateHz);
    if (sample < 0)
      return 0.0;
    if (sample >= (int)values_.size())
      return 0.0;
    return values_[sample];
  }

  // Raw sample
  float at(int i) const {
    if (i < 0 || i >= (int)values_.size())
      return 0.0f;
    return values_[i];
  }

private:
  const Point mic_pos_;
  std::vector<float> values_;
};

// Three sine waves.
static float sound_generator1(float t) {
  return 0.2 * (sin(kTestSourceFrequency * t * tau) +
                sin(kTestSourceFrequency / 3 * t * tau) +
                sin(kTestSourceFrequency / 5 * t * tau));
}

// Some sine way
static float sound_generator2(float t) {
  return 0.2 * sin(2.1637 * kTestSourceFrequency * t * tau);
}

static float sound_generator3(float t) {
  constexpr float duty_cycle = 0.1;
  return fmod(kTestSourceFrequency * (t + 1), 1) < duty_cycle ? 0.2 : -0.2;
  // return fmod(2*kTestSourceFrequency * (t + 1), 2); // sawtooth
}

uint32_t adler32(const void *buf, size_t buflength) {
  const uint8_t *buffer = (const uint8_t *)buf;

  uint32_t s1 = 1;
  uint32_t s2 = 0;

  for (size_t n = 0; n < buflength; n++) {
    s1 = (s1 + buffer[n]) % 65521;
    s2 = (s2 + s1) % 65521;
  }
  return (s2 << 16) | s1;
}

static float sound_generator4(float t) {
  uint32_t ms = (int)roundf(t * 6327);
  uint32_t v = adler32(&ms, sizeof(ms));
  return (v % 1000) / 2000.0;
}

// Correlate given microphones looking into the direction the normal
// points to.
static float cross_correlate(const MicrophoneSamples &mic1,
                             const MicrophoneSamples &mic2,
                             const Point &normal) {
  const float d1 = normal.dotMul(mic1.pos() - optical_camera_pos);
  const float d2 = normal.dotMul(mic2.pos() - optical_camera_pos);
  const float td1 = d1 / kSpeedOfSound;
  const float td2 = d2 / kSpeedOfSound;
  const int offset_1 = td1 * kSampleRateHz;
  const int offset_2 = td2 * kSampleRateHz;
  float sum = 0.0f;
  constexpr int sample_count = kMeasureTime * kSampleRateHz;
  int step = sample_count / 101;
  for (int t = 0; t < sample_count - 1000; t += step) {
    sum += mic1.at(t + offset_1) * mic2.at(t + offset_2);
  }
  return sum;
}

template <class T> class Buffer2D {
public:
  Buffer2D(int w, int h)
      : width_(w), height_(h), buffer_(new T[width_ * height_]) {
    memset(buffer_, 0, width_ * height_);
  }
  ~Buffer2D() { delete[] buffer_; }

  T &at(int x, int y) {
    assert(x >= 0 && x < width_ && y >= 0 && y < height_);
    return buffer_[y * width_ + x];
  }

  T Avg() {
    const int items = width_ * height_;
    T sum = 0;
    T *end = buffer_ + items;
    T *t = buffer_;
    while (t < end)
      sum += *t++;
    return sum / items;
  }

private:
  const int width_;
  const int height_;
  T *buffer_;
};

std::vector<Point> CreateSpiralMicrophones(int count, float base_radius,
                                           int spirals, float pitch) {
  std::vector<Point> result;
  for (int i = 0; i < count; ++i) {
    const int spiral_num = i % spirals;
    const float microphone_fraction = 1.0 * i / count;
    const float mic_radius = microphone_fraction * pitch + base_radius;
    Point microphone_pos = {
        cos(spiral_num * tau / spirals + tau / count * i / spirals) *
            mic_radius,
        sin(spiral_num * tau / spirals + tau / count * i / spirals) *
            mic_radius,
        0};
    result.push_back(microphone_pos);
  }
  return result;
}

std::vector<Point> CreateGridMicrophones(int count, float grid_dist) {
  std::vector<Point> result;
  const int x_mics = (int)sqrt(count);
  for (int i = 0; i < count; ++i) {
    int mx = i % x_mics;
    int my = i / x_mics;
    fprintf(stderr, " (%d,%d)", mx, my);
    Point microphone_pos = {mx * grid_dist, my * grid_dist, 0};
    result.push_back(microphone_pos);
  }
  return result;
}

int main(int argc, char *argv[]) {
  constexpr int sample_count = kMeasureTime * kSampleRateHz;
  const int microphones = 12;
  const float base_radius = 60 cm;
  float fmin = 1e9; // Microphone range minimum/maximum. To display.
  float fmax = -1e9;
#if 1
  std::vector<Point> micro_positions =
      CreateSpiralMicrophones(microphones, base_radius, 3, 0 cm);
#else
  std::vector<Point> micro_positions =
      CreateGridMicrophones(microphones, 25 cm);
#endif

  // Samples of what each mirophone sees
  std::vector<MicrophoneSamples *> samples;

  for (const Point &micro_pos : micro_positions) {
    if (micro_pos.x < fmin)
      fmin = micro_pos.x;
    if (micro_pos.y < fmin)
      fmin = micro_pos.y;
    if (micro_pos.x > fmax)
      fmax = micro_pos.x;
    if (micro_pos.y > fmax)
      fmax = micro_pos.y;

    MicrophoneSamples *recording =
        new MicrophoneSamples(micro_pos, sample_count);
    float distance = micro_pos.distance_to(sound_source1);
    recording->fill(sound_generator4, distance);

    distance = micro_pos.distance_to(sound_source2);
    recording->fill(sound_generator2, distance);

    samples.push_back(recording);
  }

  fprintf(stderr, "yaw: +/- %.1f [X]; pitch: +/- %.1f [Y]\n\n",
          display_range / 2 * rad2deg, display_range / 2 * rad2deg);
  TerminalCanvas canvas(resolution_pixels, resolution_pixels);

  // Print microphones.
  for (const Point &micro_pos : micro_positions) {
    int x = int(resolution_pixels * (micro_pos.x - fmin) / (fmax - fmin));
    int y = int(resolution_pixels * (micro_pos.y - fmin) / (fmax - fmin));
    canvas.SetPixel(x, y, 255, 255, 255);
  }
  canvas.Send(STDOUT_FILENO, false);

  const float range = tan(display_range / 2); // max x in one meter

  // Cross correlate
  Buffer2D<float> frame_buffer(resolution_pixels, resolution_pixels);
  long correlations = 0;
  for (int x = 0; x < resolution_pixels; ++x) {
    for (int y = 0; y < resolution_pixels; ++y) {
      // From our place, determine the vector where we're looking at.
      const float xpix = range * x / resolution_pixels - range / 2;
      const float ypix = range * y / resolution_pixels - range / 2;
      Point listen_dir = {xpix, ypix, 1};
      listen_dir.MakeUnitLen();
      float v = 0;
      for (int i = 0; i < microphones; ++i) {
        for (int j = i + 1; j < microphones; ++j) {
          v += cross_correlate(*samples[i], *samples[j], listen_dir);
          ++correlations;
        }
      }
      // The way angles are calculated from right to left, but our
      // x going from left to right, we have to mirror it.
      frame_buffer.at(resolution_pixels - x - 1, y) = v;
    }
  }

  fprintf(stderr,
          "\n%d mics; %.1f cm view in 1 meter; r=%.1fcm; f₀=%.0f; "
          "λ=%.2f cm; %.3fM cross-correlations\n",
          microphones, range * 100, base_radius * 100, kTestSourceFrequency,
          kSpeedOfSound / kTestSourceFrequency * 100, correlations / 1e6);

  float smallest = 1e9;
  float biggest = -1e9;
  for (int x = 0; x < resolution_pixels; ++x) {
    for (int y = 0; y < resolution_pixels; ++y) {
      float v = frame_buffer.at(x, y);
      if (v < smallest)
        smallest = v;
      if (v > biggest)
        biggest = v;
    }
  }

  // smallest = frame_buffer.Avg();

  const int colormap_entries = arraysize(kColorMap);
  float yaw_angle, pitch_angle;
  for (int x = 0; x < resolution_pixels; ++x) {
    for (int y = 0; y < resolution_pixels; ++y) {
      float v = frame_buffer.at(x, y);
      if (v < smallest) {
        canvas.SetPixel(x, y, 0, 0, 0);
        continue;
      }
      int color_index =
          (int)((colormap_entries - 1) * (v - smallest) / (biggest - smallest));
      if (color_index < 0) {
        fprintf(stderr, "%d\n", color_index);
        color_index = 0;
      }
      if (color_index > 255) {
        fprintf(stderr, "%d\n", color_index);
        color_index = 255;
      }
      const RGBCol &color = kColorMap[color_index];
      canvas.SetPixel(x, y, color.r * 255, color.g * 255, color.b * 255);
    }
  }

#if 1
  // it looks like, we're a little bit off. Might just be the perspective.
  Point s = sound_source1;
  s.MakeUnitLen();
  canvas.SetPixel((s.x / range + 0.5) * resolution_pixels,
                  ((0 - s.y) / range + 0.5) * resolution_pixels, 255, 255, 255);
  s = sound_source2;
  s.MakeUnitLen();
  canvas.SetPixel((s.x / range + 0.5) * resolution_pixels,
                  ((0 - s.y) / range + 0.5) * resolution_pixels, 255, 255, 255);
#endif

  canvas.Send(STDOUT_FILENO, false);
}
