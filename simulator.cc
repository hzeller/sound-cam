/*
 * Simulation of sound cam
 */

#include <assert.h>
#include <vector>
#include <functional>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <string.h>

#include "terminal-canvas.h"
#include "colormap.h"

constexpr float tau = 2 * M_PI;
constexpr float rad2deg = 360 / tau;

constexpr float kSpeedOfSound = 340.0f; // m/s
constexpr int kNoiseBits = 2;
#if 0
// Regular
constexpr int kSampleBits = 12;     
constexpr int kSampleRateHz = 44100;  // Hz
#else
// high time resolution.
constexpr int kSampleBits = 12;     
constexpr int kSampleRateHz = 192000;  // Hz
#endif

//constexpr int kSampleRateHz = 1000000;
constexpr float kMeasureTime = 1;     // second;
constexpr float kTestSourceFrequency = 1200.0;

constexpr float display_range = tau / 4;

constexpr int resolution_pixels = 101;

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

    float length() const { return sqrtf(square(x)+square(y)+square(z)); }
    
    float dotMul(const Point &other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    void MakeUnit() {
        float len = length();
        x /= len; y /= len; z /= len;
    }
    
    float x, y, z;
};

Point operator + (const Point &a, const Point &b) {
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}
Point operator - (const Point &a, const Point &b) {
    return { a.x - b.x, a.y - b.y, a.z - b.z };
}
Point operator * (const Point &p, float scalar) {
    return { p.x * scalar, p.y * scalar, p.z * scalar };
}

Point sound_source1 = { 0.3, 0.3, 1 };
Point sound_source2 = { 0, -0.4, 1 };

Point optical_camera_pos = { 0, 0, 0 };

#define arraysize(a) sizeof(a) / sizeof(a[0])

// Sound source. Returns samples from -1 to 1 at given time.
typedef std::function<float(float)> SoundSource;

class MicrophoneSamples {
public:
    MicrophoneSamples(const Point &microphone,
                      int samples) : mic_pos_(microphone), values_(samples) {}

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
            const int value = quant * source(i * sample_distance - time_delay)
                + random() % noisefloor_magnitude;
            values_[i] += (1.0 * value / quant);
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
        if (i < 0 || i >= (int)values_.size()) return 0.0f;
        return values_[i];
    }
    
private:
    const Point mic_pos_;
    std::vector<float> values_;
};

static float sound_generator1(float t) {
    return 0.1 * (sin(kTestSourceFrequency * t * tau) +
                  sin(kTestSourceFrequency/3 * t * tau) +
                  sin(kTestSourceFrequency/5 * t * tau));
}

static float sound_generator2(float t) {
    return 0.2 * sin(2.1637 * kTestSourceFrequency * t * tau);
    constexpr float duty_cycle = 0.1;
    return fmod(kTestSourceFrequency * (t + 1), 1) < duty_cycle ? 0.2 : -0.2;
    //return fmod(2*kTestSourceFrequency * (t + 1), 2); // sawtooth
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
    constexpr int sample_count = kMeasureTime * kSampleRateHz;
    const int offset_1 = td1 * kSampleRateHz;
    const int offset_2 = td2 * kSampleRateHz;    
    float sum = 0.0f;
    int step = sample_count / 100;
    for (int t = 0; t < sample_count; t += step) {
        sum += mic1.at(t + offset_1) * mic2.at(t + offset_2);
    }
    return sum;
}

template<class T> class Buffer2D {
public:
    Buffer2D(int w, int h) : width_(w), height_(h),
                             buffer_(new T [width_ * height_]) {
        memset(buffer_, 0, width_*height_);
    }
    ~Buffer2D() { delete [] buffer_; }
    
    T& at(int x, int y) {
        assert(x >= 0 && x < width_ && y >= 0 && y < height_);
        return buffer_[y * width_ + x];
    }
    
private:
    const int width_;
    const int height_;
    T *buffer_;
};

int main(int argc, char *argv[]) {
    constexpr int sample_count = kMeasureTime * kSampleRateHz;
    const int microphones = 25;
    const int spirals = 3;
    const float base_radius = 60 cm;
    float fmin = 1e9;
    float fmax = -1e9;
    std::vector<MicrophoneSamples*> samples(microphones);
    for (int i = 0; i < microphones; ++i) {
        //fprintf(stderr, "Microphone %d: ", i);
#if 1
        const int spiral_num = i % spirals;
        const float microphone_fraction = 1.0 * i / microphones;
        const float mic_radius = microphone_fraction * 0 cm + base_radius;
        //const float mic_radius = 33 cm;
        Point microphone_pos = { cos(spiral_num * tau/spirals + tau/microphones*i/spirals) * mic_radius,
                                 sin(spiral_num * tau/spirals + tau/microphones*i/spirals) * mic_radius,
                                 0 };
#else
        const int x_mics = sqrt(microphones);
        int my = i / x_mics;
        int mx = i % x_mics;
        //fprintf(stderr, " (%d,%d)", mx, my);
        Point microphone_pos = { mx * 25 cm, my * 25 cm, 0 };
        //Point microphone_pos = { rand() % 60 cm, rand() % 60 cm, 0 };
        //Point microphone_pos = { (mx * 10 + (rand() % 6 - 3)) cm, (my * 10 + (rand() % 6 - 3)) cm, 0 };
#endif
        if (microphone_pos.x < fmin) fmin = microphone_pos.x;
        if (microphone_pos.y < fmin) fmin = microphone_pos.y;
        if (microphone_pos.x > fmax) fmax = microphone_pos.x;
        if (microphone_pos.y > fmax) fmax = microphone_pos.y;
        samples[i] = new MicrophoneSamples(microphone_pos, sample_count);
        float distance = microphone_pos.distance_to(sound_source1);
        //fprintf(stderr, " %.2f cm ", distance * 100);
        samples[i]->fill(sound_generator1, distance);
        distance = microphone_pos.distance_to(sound_source2);
        samples[i]->fill(sound_generator2, distance);
        //fprintf(stderr, " %.2f cm\n", distance * 100);
    }

#if 1
    FILE *a = fopen("/tmp/a", "w");
    for (int i = 0; i < sample_count / kTestSourceFrequency; ++i) {
        fprintf(a, "%d %.4f\n", i, samples[0]->at(i));
    }
    fclose(a);
    FILE *b = fopen("/tmp/b", "w");
    for (int i = 0; i < sample_count / kTestSourceFrequency; ++i) {
        fprintf(b, "%d %.4f\n", i, samples[1]->at(i));
    }
    fclose(b);    
#endif
    
#if 1
    FILE *c = fopen("/tmp/c", "w");
    auto start = std::chrono::system_clock::now();
    float max_x = tan(display_range/2);  // max x in one meter
    for (float x = -max_x; x < max_x; x+=(2*max_x)/resolution_pixels) {
        Point listen_direction = { x, 0, 1 };
        listen_direction.MakeUnit();
        fprintf(c, "%.2f %.3f\n", x,
                cross_correlate(*samples[0], *samples[1], listen_direction));
    }
    auto duration = std::chrono::system_clock::now() - start;
    std::cout << "Duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << "ms\n";
    fclose(c);
#endif

    fprintf(stderr, "yaw: +/- %.1f [X]; pitch: +/- %.1f [Y]\n\n",
            display_range / 2 * rad2deg, display_range / 2 * rad2deg);
    TerminalCanvas canvas(resolution_pixels, resolution_pixels);

    for (int i = 0; i < microphones; ++i) {
        int x = int(resolution_pixels * (samples[i]->pos().x - fmin) / (fmax-fmin));
        int y = int(resolution_pixels * (samples[i]->pos().y - fmin) / (fmax-fmin));
        canvas.SetPixel(x, y, 255, 255, 255);
    }
    canvas.Send(STDOUT_FILENO, false);
    
    Buffer2D<float> frame_buffer(resolution_pixels, resolution_pixels);
    const float range = tan(display_range/2);  // max x in one meter
    fprintf(stderr, "\n%d mics; %.1f cm view in 1 meter; r=%.1fcm; f₀=%.0f; "
            "λ=%.2f cm\n",
            microphones, range * 100, base_radius * 100, kTestSourceFrequency,
            kSpeedOfSound / kTestSourceFrequency * 100);
    for (int x = 0; x < resolution_pixels; ++x) {
        for (int y = 0; y < resolution_pixels; ++y) {
            const float xpix = range * x / resolution_pixels - range/2;
            const float ypix = range * y / resolution_pixels - range/2;
            Point listen_dir = { xpix, ypix, 1 };
            listen_dir.MakeUnit();
            float v = 0;
            for (int i = 0; i < microphones; ++i) {
                for (int j = i+1; j < microphones; ++j) {
                    v += cross_correlate(*samples[i], *samples[j], listen_dir);
                }
            }
            // The way angles are calculated from right to left, but our
            // x going from left to right, we have to mirror it.
            frame_buffer.at(resolution_pixels - x - 1, y) = v;
        }
    }

    float smallest = 1e9;
    float biggest = -1e9;
    for (int x = 0; x < resolution_pixels; ++x) {
        for (int y = 0; y < resolution_pixels; ++y) {
            float v = frame_buffer.at(x, y);
            if (v < smallest) smallest = v;
            if (v > biggest) biggest = v;
        }
    }

    const int colormap_entries = arraysize(kColorMap);    
    float yaw_angle, pitch_angle;
    for (int x = 0; x < resolution_pixels; ++x) {
        for (int y = 0; y < resolution_pixels; ++y) {
            int color_index = (int)(
                (colormap_entries-1)
                * (frame_buffer.at(x, y) - smallest) / (biggest-smallest));
            if (color_index < 0) {
                fprintf(stderr, "%d\n", color_index);
                color_index = 0;
            }
            if (color_index > 255) {
                fprintf(stderr, "%d\n", color_index);
                color_index = 255;
            }
            const RGBCol &color = kColorMap[color_index];
            canvas.SetPixel(x, y,
                            color.r * 255, color.g * 255, color.b * 255);
        }
    }

#if 1
    // it looks like, we're a little bit off
    canvas.SetPixel((sound_source1.x/range + 0.5) * resolution_pixels,
                    ((0-sound_source1.y)/range + 0.5) * resolution_pixels,
                    255, 255, 255);

    canvas.SetPixel((sound_source2.x/range + 0.5) * resolution_pixels,
                    ((0-sound_source2.y)/range + 0.5) * resolution_pixels,
                    255, 255, 255);
#endif
    
    canvas.Send(STDOUT_FILENO, false);

    Point p1 = { -50 cm, 0, 0 };
    Point p2 = { 50 cm, 0, 0 };
    Point pc = { 0, 0, 0 };
    Point n = Point({ 5 cm, 0, 1 });
    n.MakeUnit();
    const float d1 = n.dotMul(p1 - pc);
    const float d2 = n.dotMul(p2 - pc);
    fprintf(stderr, "Distance: %.3f, %.3f  %.3f\n", d1, d2, d1-d2);
}
