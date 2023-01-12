#ifndef BUFFER_2D_H
#define BUFFER_2D_H

#include <cassert>

template <class T> class Buffer2D {
public:
  Buffer2D(int w, int h)
      : width_(w),
        height_(h),
        buffer_(new T[width_ * height_]),
        buffer_end_(buffer_ + width_ * height_) {}
  ~Buffer2D() { delete[] buffer_; }

  T &at(int x, int y) {
    assert(x >= 0 && x < width_ && y >= 0 && y < height_);
    return buffer_[y * width_ + x];
  }

  const T &at(int x, int y) const {
    assert(x >= 0 && x < width_ && y >= 0 && y < height_);
    return buffer_[y * width_ + x];
  }

  int width() const { return width_; }
  int height() const { return height_; }

  const T *begin() const { return buffer_; }
  const T *end() const { return buffer_end_; }

private:
  const int width_;
  const int height_;
  T *const buffer_;
  const T *const buffer_end_;
};

#endif // BUFFER_2D_H
