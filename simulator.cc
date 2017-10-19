/*
 * Simulation of sound cam
 */

#include <vector>
#include <functional>
#include <stdint.h>
#include <math.h>
#include <unistd.h>

#include "terminal-canvas.h"

constexpr float tau = 2 * M_PI;
constexpr float rad2deg = 360 / tau;

constexpr float kSpeedOfSound = 340.0f; // m/s
#if 0
// Regular
constexpr int kSampleBits = 15;     
constexpr int kSampleRateHz = 44100;  // Hz
#else
// Super precise
constexpr int kSampleBits = 24;     
constexpr int kSampleRateHz = 192000;  // Hz
#endif

//constexpr int kSampleRateHz = 1000000;
constexpr float kMeasureTime = 1;     // second;
constexpr float kTestSourceFrequency = 3000.0;

constexpr float display_range = tau / 4;

constexpr int resolution_pixels = 80;

#define cm * 0.01

static float square(float a) { return a * a; }
class Point {
public:
    static Point FromAngle(float yaw, float pitch) {
        return { cos(yaw) * cos(pitch),
                sin(yaw) * cos(pitch),
                sin(pitch) };
    }
        
    float distance_to(const Point& other) const {
        return sqrtf(square(x - other.x)
                     + square(y - other.y)
                     + square(z - other.z));
    }
    
    float x, y, z;
};

Point operator + (const Point &a, const Point &b) {
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}
Point operator * (const Point &p, float scalar) {
    return { p.x * scalar, p.y * scalar, p.z * scalar };
}

Point microphone_pos[4] = {
    { 0, 5 cm,  5 cm },
    { 0, -5 cm,  5 cm },
    { 0, -5 cm, -5 cm },
    { 0, 5 cm,  -5 cm },
};

//Point sound_source = { 20 cm, 0 cm, 5 cm };
Point sound_source = Point::FromAngle(10 / rad2deg, 0) * 10 cm;

Point optical_camera_pos = { 0, 0, 0 };

#define arraysize(a) sizeof(a) / sizeof(a[0])

// Sound source. Returns samples from -1 to 1 at given time.
typedef std::function<float(float)> SoundSource;

class MicrophoneSamples {
public:
    MicrophoneSamples(const Point &microphone) : mic_pos_(microphone) {}

    const Point &pos() const { return mic_pos_; }
    
    void fill(int samples, SoundSource source, float distance_m) {
        const float sample_distance = 1.0f / kSampleRateHz;
        const float time_delay = distance_m / kSpeedOfSound;
        values_.clear();
        // We keep things symmetric around 0, so half above, half below; -1bit
        const int quant = 1 << (kSampleBits - 1);
        for (int i = 0; i < samples; ++i) {
            const int value = quant * source(i * sample_distance - time_delay);
            values_.push_back(1.0 * value / quant);
        }
    }

    float value_at(float t) {
        // Simplest: find closest to sample point
        int sample = round(t * kSampleRateHz);
        if (sample < 0) return 0.0;
        if (sample >= (int)values_.size()) return 0.0;
        return values_[sample];
    }

    // Raw sample
    float at(int i) const {
        if (i < 0 || i >= values_.size()) return 0.0f;
        return values_[i];
    }
    
private:
    const Point mic_pos_;
    std::vector<float> values_;
};

static float sound_generator(float t) {
    //return sin(kTestSourceFrequency * t * tau);
    return 0.03 * (sin(kTestSourceFrequency * t * tau) +
                  sin(kTestSourceFrequency/3 * t * tau) +
                  sin(kTestSourceFrequency/5 * t * tau));
    //return fmod(2 * kTestSourceFrequency * t, 2) < 1 ? -1 : 1;
}

static float cross_correlate(const MicrophoneSamples &mic1,
                             const MicrophoneSamples &mic2,
                             float angle) {
    // We actually want the angle on plane; this is fudged. We only
    // use one microphone and multiply by two.
    const float offset_distance = 2 * sin(angle) * mic1.pos().distance_to(optical_camera_pos);
    const float time_delay = offset_distance / kSpeedOfSound;
    constexpr int sample_count = kMeasureTime * kSampleRateHz;
    const int offset_samples = time_delay * kSampleRateHz;
    float sum = 0.0f;
    for (int t = 0; t < sample_count; ++t) {
        sum += mic1.at(t) * mic2.at(t + offset_samples);
    }
    return sum;
}

int main(int argc, char *argv[]) {
    constexpr int sample_count = kMeasureTime * kSampleRateHz;
    const int microphones = arraysize(microphone_pos);
    std::vector<MicrophoneSamples*> samples(microphones);
    for (int i = 0; i < microphones; ++i) {
        samples[i] = new MicrophoneSamples(microphone_pos[i]);
        const float distance = microphone_pos[i].distance_to(sound_source);
        fprintf(stderr, "Microphone %d: %.2f cm distance\n", i, distance * 100);
        samples[i]->fill(sample_count, sound_generator, distance);
    }

#if 1
    FILE *c = fopen("/tmp/c", "w");
    for (float angle = -display_range/2; angle < display_range/2; angle+=display_range/resolution_pixels) {
        fprintf(c, "%.3f %.3f\n", angle * rad2deg,
                cross_correlate(*samples[0], *samples[1], angle));
    }
    fclose(c);
#endif
    
#if 0
    fprintf(stderr, "yaw: +/- %.1f [X]; pitch: +/- %.1f [Y]\n\n",
            display_range / 2 * rad2deg, display_range / 2 * rad2deg);
    TerminalCanvas canvas(resolution_pixels + 10, resolution_pixels);
    float yaw_angle, pitch_angle;
    for (int x = 0; x < resolution_pixels; ++x) {
        for (int y = 0; y < resolution_pixels; ++y) {
            yaw_angle = (1.0 * x / resolution_pixels - 0.5) * display_range;
            pitch_angle = (1.0 * y / resolution_pixels - 0.5) * display_range;

#if 0
            canvas.SetPixel(x, y,
                            255 * yaw_angle / display_range + 128,
                            255 * pitch_angle / display_range + 128,
                            0);
#else
            canvas.SetPixel(x, y, 127 * cos(4 * yaw_angle) + 128,
                            127 * cos(4 * pitch_angle) + 128, 0);
#endif       
            // Determine distance to each of the microphones
        }
    }

    canvas.Send(STDOUT_FILENO, false);
#endif
}
