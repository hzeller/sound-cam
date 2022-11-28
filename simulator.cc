
#include <unistd.h>
#include <functional>
#include <random>

#include "buffer-2d.h"
#include "colormap.h"
#include "cross-correlation.h"
#include "point.h"
#include "terminal-canvas.h"

// Random Microphone arrangement can reduce regular repeating artifacts
//#define USE_RANDOM_MICROPHONE_ARRANGEMENT

constexpr float tau = 2 * M_PI;
constexpr float kSpeedOfSound = 340.0f; // m/s
constexpr float kTestSourceFrequency = 1200.0; // baseline sound frequency.

constexpr int screen_size = 100;

constexpr float kMicrophoneRadius = 0.3;

constexpr Point optical_camera_pos = {0, 0, 0};

constexpr int kMicrophoneCount = 17;

constexpr size_t kSampleRateHz = 48000;
constexpr size_t kMicrophoneSamples = 1 << 9;

constexpr float display_range = tau / 4; // Angle of view. 90 degree.

typedef std::vector<float> MicrophoneRecording;
typedef std::vector<float> CrossCorrelation;
typedef std::function<float(float t)> WaveExpr;

#define arraysize(a) sizeof(a) / sizeof(a[0])

static float sampling_noise() {
  static std::random_device r_engine;
  static std::default_random_engine r(r_engine());

  // Additional noise to signal (which is in the range from -1 to 1).
  static std::uniform_real_distribution<> distribution(-1, 1);

  return distribution(r);
}

// Various microphone arrangements.
void AddMicrophoneCircle(std::vector<Point> *mics, int count, float radius) {
  fprintf(stderr, "Circle microphone arrangement\n");
  for (int i = 0; i < count; ++i) {
    const float angle = tau / count * i;
    mics->push_back(Point{cos(angle) * radius, sin(angle) * radius, 0});
  }
}

void AddMicrophoneRandom(std::vector<Point> *mics, int count, float radius) {
  fprintf(stderr, "Random microphone arrangement\n");
  std::vector<Point> microphones;
  srandom(time(NULL));
  for (int i = 0; i < count; ++i) {
    const float x = (random() % 2000 - 1000) / 1000.0 * radius;
    const float y = (random() % 2000 - 1000) / 1000.0 * radius;
    mics->push_back(Point{x, y, 0});
  }
}

// Slightly different frequencies for the wave generating functions to be
// able to distinguish them easily.
float wave1(float t) {
  return sin(2 * kTestSourceFrequency * t * tau);
};

static float wave2(float t) {
  return sin(2.1637 * kTestSourceFrequency * t * tau);
}

static float wave3(float t) {
  return sin(2.718 * kTestSourceFrequency * t * tau);
}

static const std::pair<Point, WaveExpr> sound_sources[] = {
  {{-0.8, 0.2, 4}, wave1},
  {{0, -8, 20}, wave2},
  {{1, 0, 2}, wave3},
};

// Add a recording with the given phase shift and wave.
void add_recording(MicrophoneRecording *recording,
                   int sample_frequency_hz,
                   float phase_shift_seconds,
                   std::function<float(float t)> wave_f) {
  for (size_t i = 0; i < recording->size(); ++i) {
    const float t = phase_shift_seconds + 1.0f * i / sample_frequency_hz;
    (*recording)[i] += wave_f(t) + sampling_noise();
  }
}

void VisualizeMicrophoneLocations(const std::vector<Point>& microphones) {
  float fmin = 0, fmax = 0;
  for (const Point &micro_pos : microphones) {
    if (micro_pos.x < fmin) fmin = micro_pos.x;
    if (micro_pos.y < fmin) fmin = micro_pos.y;
    if (micro_pos.x > fmax) fmax = micro_pos.x;
    if (micro_pos.y > fmax) fmax = micro_pos.y;
  }

  // Print microphones.
  fprintf(stderr, "Microphone min=%.2f max=%.2f\n", fmin, fmax);
  TerminalCanvas canvas(screen_size, screen_size);
  for (const Point &m : microphones) {
    int x = int((screen_size - 1) * (m.x - fmin) / (fmax - fmin));
    int y = int((screen_size - 1) * (m.y - fmin) / (fmax - fmin));
    canvas.SetPixel(x, y, 255, 255, 255);
  }
  canvas.Send(STDOUT_FILENO, false);
}

int main() {
  // Maximum interesting correlation time difference needed in cross-correlation
  const int kCrossCorrelateElementsOfInterest = 75;

  std::vector<Point> microphones;
#ifdef USE_RANDOM_MICROPHONE_ARRANGEMENT
  AddMicrophoneRandom(&microphones, kMicrophoneCount, kMicrophoneRadius);
  VisualizeMicrophoneLocations(microphones);
#else
  AddMicrophoneCircle(&microphones, kMicrophoneCount, kMicrophoneRadius);
#endif
  
  std::vector<MicrophoneRecording> microphone_recording;

  // Create a recording that each microphone sees.
  // That is dependent on the distance from the sound source.
  for (const Point &micro_pos : microphones) {
    MicrophoneRecording recording(kMicrophoneSamples);
    for (const auto &s : sound_sources) {
      float distance = micro_pos.distance_to(s.first);
      add_recording(&recording, kSampleRateHz, distance / kSpeedOfSound, s.second);
    }
    microphone_recording.push_back(recording);
  }

  // Prepare cross correlations between each microphone pair and remember them
  // for quick look-up later.
  Buffer2D<CrossCorrelation> microphone_cross_correlation(kMicrophoneCount,
                                                          kMicrophoneCount);
  int correlation_count = 0;
  for (int i = 0; i < kMicrophoneCount; ++i) {
    for (int j = 0; j < kMicrophoneCount; ++j) {
      if (i == j) continue;  // no need to auto-correlate.
      microphone_cross_correlation.at(i, j)
        = cross_correlate(microphone_recording[i],
                          microphone_recording[j],
                          kMicrophoneSamples/2,
                          kCrossCorrelateElementsOfInterest);
      ++correlation_count;
    }
  }
  fprintf(stderr, "Created %d cross correlations (with each %d output count)\n",
          correlation_count, kCrossCorrelateElementsOfInterest);

  // Sweep.
  const float range = std::tan(display_range / 2); // max x in one meter
  Buffer2D<float> frame_buffer(screen_size, screen_size);
  int max_offset_used = 0;
  for (int x = 0; x < screen_size; ++x) {
    for (int y = 0; y < screen_size; ++y) {
      // From our place, determine the vector where we're looking at.
      const float xpix = range * x / screen_size - range / 2;
      const float ypix = range * y / screen_size - range / 2;
      Point listen_dir = {xpix, ypix, 1};
      listen_dir.MakeUnitLen();
      float value = 0;
      // Determine distance of plane facing where we're looking at from
      // microphone. Do that for each pair of microphones and look-up the
      // corresponding cross correlation.
      for (int i = 0; i < kMicrophoneCount; ++i) {
        const float d1 = listen_dir.dotMul(microphones[i] - optical_camera_pos);
        const float td1 = d1 / kSpeedOfSound;
        for (int j = i + 1; j < kMicrophoneCount; ++j) {
          const float d2 = listen_dir.dotMul(microphones[j] - optical_camera_pos);
          const float td2 = d2 / kSpeedOfSound;
          const int offset = roundf((td2 - td1) * kSampleRateHz);
          assert(abs(offset) < kCrossCorrelateElementsOfInterest);
          if (offset >= 0) {
            value += microphone_cross_correlation.at(i, j)[offset];
          }
          else {
            value += microphone_cross_correlation.at(j, i)[-offset];
          }
          if (abs(offset) > max_offset_used) max_offset_used = abs(offset);
        }
      }
      // The way angles are calculated from right to left, but our
      // x going from left to right, we have to mirror it.
      frame_buffer.at(screen_size - x - 1, y) = value;
    }
  }
  fprintf(stderr, "Maximum cross-correlate output count used: %d\n",
          max_offset_used);

  fprintf(stderr,
          "\n%d mics; %.1f cm view in 1 meter; r=%.1fcm; f₀=%.0f; "
          "λ=%.2f cm; %.3fms max offset\n",
          kMicrophoneCount, range * 100, kMicrophoneRadius * 100,
          kTestSourceFrequency,
          kSpeedOfSound / kTestSourceFrequency * 100,
          max_offset_used * 1000.0 / kSampleRateHz);

  // Determine range for the coloring.
  float smallest = 1e9;
  float biggest = -1e9;
  for (int x = 0; x < screen_size; ++x) {
    for (int y = 0; y < screen_size; ++y) {
      float v = frame_buffer.at(x, y);
      if (v < smallest)
        smallest = v;
      if (v > biggest)
        biggest = v;
    }
  }

  TerminalCanvas canvas(screen_size, screen_size);
  const int colormap_entries = arraysize(kColorMap);
  for (int x = 0; x < screen_size; ++x) {
    for (int y = 0; y < screen_size; ++y) {
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

  // Some overlay where the sound sources are.
  for (const auto &s : sound_sources) {
    Point loc = s.first;
    loc.MakeUnitLen();
    canvas.SetPixel((loc.x / range + 0.5) * screen_size,
                    ((0 - loc.y) / range + 0.5) * screen_size, 255, 255, 255);
  }

  canvas.Send(STDOUT_FILENO, false);
}
