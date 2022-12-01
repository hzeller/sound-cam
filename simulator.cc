
#include <signal.h>
#include <unistd.h>
#include <sys/time.h>
#include <termios.h>
#include <sys/select.h>
#include <fcntl.h>

#include <ctime>
#include <functional>
#include <random>

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

typedef std::vector<real_t> MicrophoneRecording;
typedef std::vector<real_t> CrossCorrelation;
typedef std::function<real_t(real_t t)> WaveExpr;

#define arraysize(a) sizeof(a) / sizeof(a[0])

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
real_t wave1(real_t t) { return sin(2 * kTestSourceFrequency * t * tau); };

static real_t wave2(real_t t) {
  return sin(2.1637 * kTestSourceFrequency * t * tau);
}

static real_t wave3(real_t t) {
  return sin(2.718 * kTestSourceFrequency * t * tau);
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

void VisualizeMicrophoneLocations(const std::vector<Point> &microphones) {
  real_t fmin = 0, fmax = 0;
  for (const Point &micro_pos : microphones) {
    if (micro_pos.x < fmin)
      fmin = micro_pos.x;
    if (micro_pos.y < fmin)
      fmin = micro_pos.y;
    if (micro_pos.x > fmax)
      fmax = micro_pos.x;
    if (micro_pos.y > fmax)
      fmax = micro_pos.y;
  }

  // Print microphones.
  fprintf(stderr, "Microphone min=%.2f max=%.2f\n", fmin, fmax);
  TerminalCanvas canvas(kScreenSize, kScreenSize);
  for (const Point &m : microphones) {
    int x = int((canvas.width() - 1) * (m.x - fmin) / (fmax - fmin));
    int y = int((canvas.height() - 1) * (m.y - fmin) / (fmax - fmin));
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
std::vector<correlate_preprocessed_t>
SimulateRecording(const std::vector<Point> &microphone_locations) {
  // That is dependent on the distance from the sound sources and is a sum
  // of these.
  std::vector<correlate_preprocessed_t> result;
  for (const Point &micro_pos : microphone_locations) {
    MicrophoneRecording recording(kMicrophoneSamples);
    for (const auto &s : sound_sources) {
      const real_t distance = micro_pos.distance_to(s.loc);
      add_recording(&recording, kSampleRateHz, distance / kSpeedOfSound, s.gen);
    }
    result.push_back(PreprocessCorrelate(recording));
  }
  return result;
}

// Prepare cross correlations between each microphone pair and remember them
// for quick look-up later.
// TODO: we might only need to one diagonal half if we swap the locations
// at look-up.
void PrecalculateCrossCorrelationMatrix(
    const std::vector<correlate_preprocessed_t>& microphone_recording,
    Buffer2D<CrossCorrelation> *cross_correlations) {
  assert((int)microphone_recording.size() == cross_correlations->width());
  assert((int)microphone_recording.size() == cross_correlations->height());
  int correlation_count = 0;
  for (size_t i = 0; i < microphone_recording.size(); ++i) {
    for (size_t j = 0; j < microphone_recording.size(); ++j) {
      if (i == j)
        continue; // no need to self cross-correlate.
      cross_correlations->at(i, j) = cross_correlate(
        microphone_recording[i], microphone_recording[j]);
      ++correlation_count;
    }
  }
#if 0
  fprintf(stderr, "Created %d cross correlations\n",
          correlation_count);
#endif
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

// Construct sound image. Sweeping a vector 'looking' into a particular
// direction, determining what the expected time difference is for each
// microphone-pair and adding up the corresponding cross correlations for
// each pixel.
void ConstructSoundImage(
    const Point &view_origin, real_t range, //
    const std::vector<Point> &microphone_loc,
    const Buffer2D<CrossCorrelation> &microphone_cross_correlation,
    Buffer2D<real_t> *frame_buffer) {

  const int microphone_count = microphone_loc.size();
  int max_offset_used = 0;
  for (int x = 0; x < frame_buffer->width(); ++x) {
    for (int y = 0; y < frame_buffer->height(); ++y) {
      // From our place, determine the vector where we're looking at.
      const real_t xpix = range * x / frame_buffer->width() - range / 2;
      const real_t ypix = range * y / frame_buffer->height() - range / 2;
      Point listen_dir = {xpix, ypix, 1};
      listen_dir.MakeUnitLen();  // normal vector of wavefront plane.
      real_t value = 0;
      // Determine distance of plane facing where we're looking at from
      // microphone. Do that for each pair of microphones and look-up the
      // corresponding cross correlation.
      for (int i = 0; i < microphone_count; ++i) {
        const real_t d1 = listen_dir.dotMul(microphone_loc[i] - view_origin);
        const real_t td1 = d1 / kSpeedOfSound;

        for (int j = i + 1; j < microphone_count; ++j) {
          const real_t d2 = listen_dir.dotMul(microphone_loc[j] - view_origin);
          const real_t td2 = d2 / kSpeedOfSound;
          const int offset = (td2 - td1) * kSampleRateHz;
          if (offset >= 0) {
            value += microphone_cross_correlation.at(i, j)[offset];
          } else {
            value += microphone_cross_correlation.at(j, i)[-offset];
          }
          if (abs(offset) > max_offset_used)
            max_offset_used = abs(offset);
        }
      }
      // The way angles are calculated from right to left, but our
      // x going from left to right, we have to mirror it.
      frame_buffer->at(frame_buffer->width() - x - 1, y) = value;
    }
  }
#if 0
  fprintf(stderr, "Maximum cross-correlate output count used: %d\n",
          max_offset_used);

  fprintf(stderr,
          "\n%d mics; %.1f cm view in 1 meter; r=%.1fcm; f₀=%.0f; "
          "λ=%.2f cm; %.3fms max offset\n",
          microphone_count, range * 100, kMicrophoneRadius * 100,
          kTestSourceFrequency, kSpeedOfSound / kTestSourceFrequency * 100,
          max_offset_used * 1000.0 / kSampleRateHz);
#endif
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

int main() {
  const auto microphones = CreateMicrophoneLocations(kMicrophoneCount);
  fprintf(stderr, "Got %d microphones\n", (int)microphones.size());

  Buffer2D<CrossCorrelation> cross_correlations(microphones.size(),
                                                microphones.size());

  Buffer2D<real_t> frame_buffer(kScreenSize, kScreenSize);

  printf("\n"
         "Highlighted source movable     |   K         |"
         "    m : show microphones\n");
  printf("1, 2, 3: choose source to move | H   L  Move |"
         " <ESC>: exit\n");
  printf("                               |   J         |"
         "    o : switch output\n");
  term_raw();

  TerminalCanvas canvas(frame_buffer.width(), frame_buffer.height());
  canvas.CursorOff(STDOUT_FILENO);

  int move_source = 0;
  bool canvas_needs_jump_to_top = false;
  size_t frame_count = 0;
  bool finished = false;
  bool do_output = true;
  const auto start_time = GetTimeInMillis();
  while (!finished) {
    // Simulate recording, including noise.
    auto microphone_recordings = SimulateRecording(microphones);

    PrecalculateCrossCorrelationMatrix(microphone_recordings,
                                       &cross_correlations);

    // Now the actual image construction
    const real_t range = std::tan(display_range / 2); // max x in one meter
    ConstructSoundImage(optical_camera_pos, range, microphones,
                        cross_correlations, &frame_buffer);

    VisualizeBuffer(frame_buffer, &canvas);
    VisualizeSoundSourceLocations(range, move_source, &canvas);

    if (do_output) {
      canvas.Send(STDOUT_FILENO, canvas_needs_jump_to_top);
    }
    canvas_needs_jump_to_top = true;
    ++frame_count;

    switch (maybe_readchar()) {
    case '\033':
      finished = true;
      break;
    case '1': move_source = 0; break;
    case '2': move_source = 1; break;
    case '3': move_source = 2; break;
    case 'h': case 'H':
      move_limited(-0.1, -1, 1, &sound_sources[move_source].loc.x);
      break;
    case 'j': case 'J':
      move_limited(-0.1, -1, 1, &sound_sources[move_source].loc.y);
      break;
    case 'l': case 'L':
      move_limited(+0.1, -1, 1, &sound_sources[move_source].loc.x);
      break;
    case 'k': case 'K':
      move_limited(+0.1, -1, 1, &sound_sources[move_source].loc.y);
      break;
    case 'm':
      VisualizeMicrophoneLocations(microphones);
      canvas_needs_jump_to_top = false;
      break;
    case 'o':
      do_output = !do_output;
      if (!do_output) fprintf(stderr, "Suspended output\n");
      break;
    }
  }
  const auto duration_ms = GetTimeInMillis() - start_time;
  fprintf(stderr, "\n%.1ffps\n", 1000.0 * frame_count / duration_ms);
  canvas.CursorOn(STDOUT_FILENO);
}
