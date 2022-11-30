#ifndef BUFFER_2D_H
#define BUFFER_2D_H

#include <cassert>

template <class T> class Buffer2D {
public:
  Buffer2D(int w, int h)
      : width_(w), height_(h), buffer_(new T[width_ * height_]) {}
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

private:
  const int width_;
  const int height_;
  T *buffer_;
};

#endif // BUFFER_2D_H
