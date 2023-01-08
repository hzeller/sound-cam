
#include <assert.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/select.h>
#include <sys/time.h>
#include <termios.h>
#include <unistd.h>

#include <ctime>
#include <functional>
#include <random>
#include <span>

#include "buffer-2d.h"
#include "colormap.h"
#include "cross-correlation.h"
#include "point.h"
#include "terminal-canvas.h"

// Random Microphone arrangement can reduce regular repeating artifacts
//#define USE_RANDOM_MICROPHONE_ARRANGEMENT

constexpr real_t tau = 2 * M_PI;
constexpr real_t kSpeedOfSound = 340.0f;        // m/s
constexpr real_t kTestSourceFrequency = 1200.0; // baseline sound frequency.

constexpr real_t kMicrophoneRadius = 0.3;

constexpr Point optical_camera_pos = {0, 0, 0};

constexpr int kScreenSize = 80; // Pixels we put on terminal screen.
constexpr int kMicrophoneCount = 17;

constexpr size_t kSampleRateHz = 48000;
constexpr size_t kMicrophoneSamples = 1 << 9;

constexpr real_t display_range = tau / 4; // Angle of view. 90 degree.

typedef complex_span_t MicrophoneRecording;
typedef std::function<real_t(real_t t)> WaveExpr;

#define arraysize(a) sizeof(a) / sizeof(a[0])

static int RoundToNextPowerOf2(int val) {
  if ((val & (val - 1)) == 0) return val;
  int bit_count = 0;
  while (val >>= 1) ++bit_count;
  return 1 << (bit_count+1);
}

// Correlation context for a pair of microphones. It operates on the already
// fixed set of input ffts to generate the cross correlation.
class PairCrossCorrelation {
public:
  // The 'constructor' (which we can't have as constructor as we initialize
  // it in a loop).
  // Params are the already pre-allocated ranges for the microphone
  // and other microphone fft.
  void InitializeMicrophonePair(const complex_span_t microphone_fft_in,
                                const complex_span_t pattern_fft_in) {
    assert(microphone_fft_in.size() == pattern_fft_in.size());
    microphone_fft = microphone_fft_in;
    pattern_fft    = pattern_fft_in;
    pair_wise_multplication_backing_store.resize(microphone_fft.size());
    cross_correlation_output_backing_store.resize(microphone_fft.size());

    pair_wise_multplication = pair_wise_multplication_backing_store;
    cross_correlation_output = cross_correlation_output_backing_store;
  }

  // Compute cross correlation.
  // Precondition: both the inputs, microphone fft, and pattern fft are already
  // computed.
  void ComputeCrossCorrelation() {
    for (size_t i = 0; i < microphone_fft.size(); ++i) {
      pair_wise_multplication[i] = microphone_fft[i] * pattern_fft[i];
    }
    InvFFT(pair_wise_multplication, &cross_correlation_output);
  }

  // Get cross correlation at given index. The reference returned will be
  // stable, i.e. after a new ComputeCrossCorrelation(), the location pointing
  // the same index will be the same (Important for the PreprocessSoundImage).
  //
  // Precondition for valid data: ComputeCrossCorrelation() has been called.
  const Complex &at(int index) const { return cross_correlation_output[index]; }

private:
  complex_span_t microphone_fft;
  complex_span_t pattern_fft;

  complex_span_t pair_wise_multplication;  // intermediate
  complex_span_t cross_correlation_output;

  // To refactor: this will come from a global storage at some point to be ready for
  // plan-many.
  complex_vec_t pair_wise_multplication_backing_store;
  complex_vec_t cross_correlation_output_backing_store;
};

class Microphone {
public:
  Microphone() {}
  Microphone(const Microphone &) = delete;
  Microphone(Microphone &&)      = delete;

  Point loc;                      // Place of the microphone
  MicrophoneRecording recording;  // samples, to be filled from source.

  complex_span_t padded_recording;  // samples with padding.
  complex_span_t microphone_fft;    // local microphone fft

  complex_span_t reverse_signal;  // Reverse samples are stored at end.
  complex_span_t pattern_fft;     // reverse signal fft

  std::vector<PairCrossCorrelation> correlation;  // correlation with the i'th other mic

  void CreateReversePatternData() {
    assert(padded_recording.size() == reverse_signal.size());
    // filling the end with reverse pattern.
    std::copy(recording.rbegin(), recording.rend(),
              reverse_signal.end() - recording.size());
  }

  void ExecuteFFTs() {
    FFT(padded_recording, &microphone_fft);
    FFT(reverse_signal, &pattern_fft);
  }

  void ClearSamples() { std::fill(recording.begin(), recording.end(), 0); }
};

typedef std::vector<Microphone> MicrophoneArray;
struct MicrophoneContainer {
  MicrophoneContainer(const std::vector<Point> &locations, int samples)
    : microphones(locations.size()),
      convolution_width(RoundToNextPowerOf2(samples * 2)),
      recording_store(microphones.size() * convolution_width),  // bulk store

      // Global allocation for intermediate and result data. These need to be
      // contiguous so that we can use it with the multi-plan fft
      reverse_store(microphones.size() * convolution_width),
      microphone_fft_store(microphones.size() * convolution_width),
      pattern_fft_store(microphones.size() * convolution_width),

      // Spacing with enough padding within the input array.
      pad_offset((convolution_width - samples) / 2) {

    fprintf(stderr, "FFT size: %d\n", (int)convolution_width);
    for (size_t i = 0; i < microphones.size(); ++i) {
      microphones[i].loc       = locations[i];
      // We have one big allocation for all microphones, but give each
      // microphone its own std::span of 'their' data. Narrower view, as it is only
      // samples wide, not the full convolution width that contains padding.
      microphones[i].recording = std::span(
        recording_store.begin() + i * convolution_width + pad_offset,
        samples);
      // ... and the span that contains the empty padding around.
      microphones[i].padded_recording = std::span(
        recording_store.begin() + i * convolution_width, convolution_width);

      // All the intermediate data pointing to a global contiguous array.
      microphones[i].reverse_signal = std::span(
        reverse_store.begin() + i * convolution_width, convolution_width);

      microphones[i].microphone_fft = std::span(
        microphone_fft_store.begin() + i * convolution_width, convolution_width);

      microphones[i].pattern_fft = std::span(
        pattern_fft_store.begin() + i * convolution_width, convolution_width);
    }

    // Initialize the pair-wise cross correlations between microphones.
    for (size_t i = 0; i < microphones.size(); ++i) {
      // In the following, we actually only need half that as we fill the
      // triangle. But for convenience we just allocate all but only fill
      // half.
      microphones[i].correlation.resize(microphones.size());
      for (size_t j = i + 1; j < microphones.size(); ++j) {
        microphones[i].correlation[j].InitializeMicrophonePair(
          microphones[i].microphone_fft, microphones[j].pattern_fft);
      }
    }
  }

  int microphone_count() const { return microphones.size(); }

  void PrepareCrossCorrelations() {
    for (Microphone &m : microphones) {
      m.CreateReversePatternData();
      m.ExecuteFFTs();   // These will eventually be a part of a multi-plan.
    }
    for (size_t i = 0; i < microphones.size(); ++i) {
      for (size_t j = i+1; j < microphones.size(); ++j) {
        microphones[i].correlation[j].ComputeCrossCorrelation();
      }
    }
  }

  // Get correlation between microphone "m1" and "m2" at sampling time offset
  const Complex &getCorrelation(size_t m1, size_t m2, int offset) const {
    constexpr int kMagicLookupOffset = -1; // unclear, why always one left ?
    if (m1 > m2) {
      std::swap(m1, m2);
      offset = -offset;
    }
    return microphones[m1]
      .correlation[m2].at(pad_offset-offset+kMagicLookupOffset);
  }

  MicrophoneArray microphones;
  const size_t convolution_width;

  // Storage of all samples of all microphones back to back with padding.
  complex_vec_t recording_store;
  complex_vec_t reverse_store;
  complex_vec_t microphone_fft_store;
  complex_vec_t pattern_fft_store;

  const int pad_offset;
};

typedef int64_t tmillis_t;
static tmillis_t GetTimeInMillis() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return tp.tv_sec * 1000 + tp.tv_usec / 1000;
}

static real_t sampling_noise() {
  static std::random_device r_engine;
  static std::default_random_engine r(r_engine());

  // Additional noise to signal (which is in the range from -1 to 1).
  static std::uniform_real_distribution<> distribution(-1.5, 1.5);

  return distribution(r);
}

// Various microphone arrangements.
void AddMicrophoneCircle(std::vector<Point> *mics, int count, real_t radius) {
  fprintf(stderr, "Circle microphone arrangement, radius %.2fm; %d mics\n",
          radius, count);
  for (int i = 0; i < count; ++i) {
    const real_t angle = tau / count * i;
    mics->push_back(Point{cos(angle) * radius, sin(angle) * radius, 0});
  }
}

void AddMicrophoneRandom(std::vector<Point> *mics, int count, real_t radius) {
  fprintf(stderr, "Random microphone arrangement\n");
  std::vector<Point> microphones;
  srandom(time(NULL));
  for (int i = 0; i < count; ++i) {
    const real_t x = (random() % 2000 - 1000) / 1000.0 * radius;
    const real_t y = (random() % 2000 - 1000) / 1000.0 * radius;
    mics->push_back(Point{x, y, 0});
  }
}

// Slightly different frequencies for the wave generating functions to be
// able to distinguish them easily and not creating cross talk.
real_t wave1(real_t t) { return sinf(2 * kTestSourceFrequency * t * tau); };

static real_t wave2(real_t t) {
  return sinf(2.1637 * kTestSourceFrequency * t * tau);
}

static real_t wave3(real_t t) {
  return sinf(2.718 * kTestSourceFrequency * t * tau);
}

// Initial placement of sound sources, but read/write as we allow to
// edit them.
static struct SoundSource {
  Point loc;
  WaveExpr gen;
} sound_sources[] = {
  {{0,       0, 1.4}, wave1},
  {{-0.2, -0.3, 1.4}, wave2},
  {{0.7,   0.3, 1.4}, wave3},
};

// Add a recording with the given phase shift and wave.
void add_recording(MicrophoneRecording *recording, int sample_frequency_hz,
                   real_t phase_shift_seconds,
                   std::function<real_t(real_t t)> wave_f) {
  for (size_t i = 0; i < recording->size(); ++i) {
    const real_t t = phase_shift_seconds + 1.0f * i / sample_frequency_hz;
    (*recording)[i] += wave_f(t) + sampling_noise();
  }
}

void VisualizeMicrophoneLocations(const MicrophoneArray &microphones) {
  real_t fmin = 0, fmax = 0;
  for (const Microphone &micro : microphones) {
    if (micro.loc.x < fmin)
      fmin = micro.loc.x;
    if (micro.loc.y < fmin)
      fmin = micro.loc.y;
    if (micro.loc.x > fmax)
      fmax = micro.loc.x;
    if (micro.loc.y > fmax)
      fmax = micro.loc.y;
  }

  // Print microphones.
  fprintf(stderr, "Microphone min=%.2f max=%.2f\n", fmin, fmax);
  TerminalCanvas canvas(kScreenSize, kScreenSize);
  for (const Microphone &m : microphones) {
    int x = int((canvas.width() - 1) * (m.loc.x - fmin) / (fmax - fmin));
    int y = int((canvas.height() - 1) * (m.loc.y - fmin) / (fmax - fmin));
    canvas.SetPixel(x, y, 255, 255, 255);
  }
  canvas.Send(STDOUT_FILENO, false);
}

std::vector<Point> CreateMicrophoneLocations(int count) {
  std::vector<Point> result;
#ifdef USE_RANDOM_MICROPHONE_ARRANGEMENT
  AddMicrophoneRandom(&result, count, kMicrophoneRadius);
#else
#if 1
  AddMicrophoneCircle(&result, count, kMicrophoneRadius);
#else
  std::initializer_list<real_t> radiuses = {0.5, 1.0, 1.3};
  real_t sum = std::accumulate(radiuses.begin(), radiuses.end(), 0);
  for (const real_t radius : radiuses) {
    AddMicrophoneCircle(&result, count * radius / sum, radius * kMicrophoneRadius);
  }
#endif
#endif
  return result;
}

// Simulate what each microphone sees.
void SimulateRecording(MicrophoneArray *microphones) {
  // That is dependent on the distance from the sound sources and is a sum
  // of these.
  for (Microphone &microphone : *microphones) {
    microphone.ClearSamples();
    for (const auto &s : sound_sources) {
      const real_t distance = microphone.loc.distance_to(s.loc);
      add_recording(&microphone.recording, kSampleRateHz,
                    distance / kSpeedOfSound, s.gen);
    }
  }
}

void VisualizeSoundSourceLocations(real_t frame_width_meter,
                                   int hightlight,
                                   TerminalCanvas *canvas) {
  // Some overlay where the sound sources are.
  for (const auto &s : sound_sources) {
    Point loc = s.loc;
    loc.MakeUnitLen();
    const int xpos = (loc.x / frame_width_meter + 0.5) * canvas->width();
    const int ypos = (-loc.y / frame_width_meter + 0.5) * canvas->height();
    if (hightlight-- == 0)
      canvas->SetPixel(xpos, ypos, 255, 255, 255);
    else
      canvas->SetPixel(xpos, ypos, 127, 127, 127);
  }
}

// Given a 2D buffer, display the available ranges of values color-coded
void VisualizeBuffer(const Buffer2D<real_t> &frame_buffer,
                     TerminalCanvas *canvas) {
  // Determine range for the coloring.
  real_t smallest = 1e9;
  real_t biggest = -1e9;
  for (int x = 0; x < canvas->width(); ++x) {
    for (int y = 0; y < canvas->height(); ++y) {
      real_t v = frame_buffer.at(x, y);
      if (v < smallest)
        smallest = v;
      if (v > biggest)
        biggest = v;
    }
  }

  const int colormap_entries = arraysize(kColorMap);
  for (int x = 0; x < canvas->width(); ++x) {
    for (int y = 0; y < canvas->height(); ++y) {
      const real_t v = frame_buffer.at(x, y);
      if (v < smallest) {
        canvas->SetPixel(x, y, 0, 0, 0);
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
      canvas->SetPixel(x, y, color.r * 255, color.g * 255, color.b * 255);
    }
  }
}

// Preprocess all the places we need to add up per pixel.
// Sweeping a vector 'looking' into a particular
// direction, determining what the expected time difference is for each
// microphone-pair and remembering the corresponding cross correlations for
// each pixel.
typedef Buffer2D<std::vector<const Complex*>> preprocess_offsets_t;
void PreprocessSoundImage(const Point &view_origin, real_t range,
                          int width, int height,
                          const MicrophoneContainer &sensor,
                          preprocess_offsets_t *addition_sites) {
  const MicrophoneArray &microphones = sensor.microphones;
  const int microphone_count = microphones.size();
  int min_offset_used = 0;
  int max_offset_used = 0;
  for (int x = 0; x < width; ++x) {
    for (int y = 0; y < height; ++y) {
      // From our place, determine the vector where we're looking at.
      const real_t xpix = range * x / width - range / 2;
      const real_t ypix = range * y / height - range / 2;
      Point listen_dir = {xpix, ypix, 1};
      listen_dir.MakeUnitLen();  // normal vector of wavefront plane.
      real_t value = 0;
      // Determine distance of plane facing where we're looking at from
      // microphone. Do that for each pair of microphones and look-up the
      // corresponding cross correlation.
      for (int i = 0; i < microphone_count; ++i) {
        const real_t d1 = listen_dir.dotMul(microphones[i].loc - view_origin);
        const real_t td1 = d1 / kSpeedOfSound;

        for (int j = i + 1; j < microphone_count; ++j) {
          const real_t d2 = listen_dir.dotMul(microphones[j].loc - view_origin);
          const real_t td2 = d2 / kSpeedOfSound;
          const int offset = (td2 - td1) * kSampleRateHz;
          addition_sites->at(x, y).push_back(&sensor.getCorrelation(j, i, offset));
          if (offset > max_offset_used)
            max_offset_used = offset;
          if (offset < min_offset_used)
            min_offset_used = offset;
        }
      }
    }
  }
  fprintf(stderr, "Maximum cross-correlate offset: %d..%d\n",
          min_offset_used, max_offset_used);
  fprintf(stderr,
          "\n%d mics; %.1f cm view in 1 meter; r=%.1fcm; f₀=%.0f; "
          "λ=%.2f cm; %.3fms max offset\n",
          microphone_count, range * 100, kMicrophoneRadius * 100,
          kTestSourceFrequency, kSpeedOfSound / kTestSourceFrequency * 100,
          max_offset_used * 1000.0 / kSampleRateHz);
}

// Construct sound image. Sweeping a vector 'looking' into a particular
// direction, determining what the expected time difference is for each
// microphone-pair and adding up the corresponding cross correlations for
// each pixel.
void ConstructSoundImage(const preprocess_offsets_t &offsets,
                         Buffer2D<real_t> *frame_buffer) {
  for (int x = 0; x < frame_buffer->width(); ++x) {
    for (int y = 0; y < frame_buffer->height(); ++y) {
      real_t value = 0;
      for (const auto &cross_correlation : offsets.at(x, y)) {
        value += cross_correlation->real();
      }
      // The way angles are calculated from right to left, but our
      // x going from left to right, we have to mirror it.
      frame_buffer->at(frame_buffer->width() - x - 1, y) = value;
    }
  }
}

static struct termios orig_term;
void term_reset() {
  fcntl(STDIN_FILENO, 0, FNDELAY);
  tcsetattr(STDIN_FILENO, TCSANOW, &orig_term);
}

void term_raw() {
  tcgetattr(STDIN_FILENO, &orig_term);
  struct termios local;
  local = orig_term;
  local.c_lflag &= ~(ECHO | ICANON);
  atexit(term_reset);
  tcsetattr(STDIN_FILENO, TCSANOW, &local);
}

// Read char if available, otherwise 0.
char maybe_readchar() {
  struct timeval tv = {0, 0};  // No wait.
  fd_set fds;
  FD_ZERO(&fds);
  FD_SET(STDIN_FILENO, &fds);
  select(STDIN_FILENO+1, &fds, nullptr, nullptr, &tv);
  if (!FD_ISSET(0, &fds))
    return 0;
  char c;
  return read(STDIN_FILENO, &c, 1) <= 0 ? 0 : c;
}

void move_limited(real_t diff, real_t min, real_t max, real_t *target) {
  real_t new_value = *target + diff;
  if (new_value > min && new_value < max)
    *target = new_value;
}

int main(int argc, char *argv[]) {
  bool simulate_recording_for_each_image = true;
  bool construct_sound_image = true;
  bool do_output = true;

  bool read_keyboard = true;
  int frame_limit = -1;

  int opt;
  while ((opt = getopt(argc, argv, "f:sr")) != -1) {
    switch (opt) {
    case 'f':
      frame_limit = atoi(optarg);
      do_output = false;  // Used to performance test, so no output.
      read_keyboard = false;
      break;
    case 's':
      construct_sound_image = false;
      break;
    case 'r':
      simulate_recording_for_each_image = false;
      break;
    default:
      fprintf(stderr, "Usage: %s [-f <frames>] [-s]\n", argv[0]);
      fprintf(stderr, " -f <frames>  : perf test: only calculate no output\n");
      fprintf(stderr, "Reduce other calculations to just focus on FFT\n");
      fprintf(stderr, " -s           : Only do FFTs, don't do calc sound\n");
      fprintf(stderr, " -r           : Only create simulated recording once\n");
      return 0;
    }
  }

  // After we have the microphone locations, we can create a pre-allocated
  // fixed set of microphones.

  MicrophoneContainer sensor(CreateMicrophoneLocations(kMicrophoneCount),
                             kMicrophoneSamples);

  fprintf(stderr, "Got %d microphones\n", sensor.microphone_count());

  Buffer2D<real_t> frame_buffer(kScreenSize, kScreenSize);

  if (read_keyboard) {
    printf(R"(
  Highlighted source movable     |   K         |    m : show microphones
  1, 2, 3: choose source to move | H   L  Move | <ESC>: exit
                                 |   J         |    o : switch output
)");
    term_raw();
  }

  TerminalCanvas canvas(frame_buffer.width(), frame_buffer.height());
  canvas.CursorOff(STDOUT_FILENO);

  const real_t range = std::tan(display_range / 2); // max x in one meter
  preprocess_offsets_t preprocessed_offsets(kScreenSize, kScreenSize);
  PreprocessSoundImage(optical_camera_pos, range,
                       frame_buffer.width(), frame_buffer.height(),
                       sensor, &preprocessed_offsets);

  int move_source = 0;
  bool canvas_needs_jump_to_top = false;
  size_t frame_count = 0;
  bool finished = false;
  bool do_image_once = true;
  const auto start_time = GetTimeInMillis();
  while (!finished && frame_limit != 0) {
    if (frame_limit > 0) --frame_limit;
    // Simulate recording, including noise.
    if (do_image_once ||
        simulate_recording_for_each_image || frame_count == 0) {
      SimulateRecording(&sensor.microphones);
    }
    sensor.PrepareCrossCorrelations();

    // Now the actual image construction
    if (construct_sound_image || do_image_once) {
      ConstructSoundImage(preprocessed_offsets, &frame_buffer);
      VisualizeBuffer(frame_buffer, &canvas);
      VisualizeSoundSourceLocations(range, move_source, &canvas);
    }

    if (do_output && (construct_sound_image || do_image_once)) {
      canvas.Send(STDOUT_FILENO, canvas_needs_jump_to_top);
    }
    canvas_needs_jump_to_top = true;
    ++frame_count;
    do_image_once = false;

    if (read_keyboard) {
        switch (maybe_readchar()) {
        case '\033': finished = true; break;
        case '1': move_source = 0; break;
        case '2': move_source = 1; break;
        case '3': move_source = 2; break;
        case 'h':
        case 'H':
            move_limited(-0.1, -1, 1, &sound_sources[move_source].loc.x);
            do_image_once = true;
            break;
        case 'j':
        case 'J':
            move_limited(-0.1, -1, 1, &sound_sources[move_source].loc.y);
            do_image_once = true;
            break;
        case 'l':
        case 'L':
            move_limited(+0.1, -1, 1, &sound_sources[move_source].loc.x);
            do_image_once = true;
            break;
        case 'k':
        case 'K':
            move_limited(+0.1, -1, 1, &sound_sources[move_source].loc.y);
            do_image_once = true;
            break;
        case 'm':
            VisualizeMicrophoneLocations(sensor.microphones);
            canvas_needs_jump_to_top = false;
            break;
        case 'o':
            do_output = !do_output;
            if (!do_output) fprintf(stderr, "Suspended output\n");
            break;
        }
    }
  }
  const auto duration_ms = GetTimeInMillis() - start_time;
  fprintf(stderr, "\n%.2ffps\n", 1000.0 * frame_count / duration_ms);
  canvas.CursorOn(STDOUT_FILENO);
}
