#ifndef BUFFER_2D_H
#define BUFFER_2D_H

#include <cassert>

template <class T> class Buffer2D {
public:
  Buffer2D(int w, int h)
    : width_(w), height_(h), buffer_(new T[width_ * height_]) {
    //memset(buffer_, 0, width_ * height_);
  }
  ~Buffer2D() { delete[] buffer_; }

  T &at(int x, int y) {
    assert(x >= 0 && x < width_ && y >= 0 && y < height_);
    return buffer_[y * width_ + x];
  }

private:
  const int width_;
  const int height_;
  T *buffer_;
};

#endif // BUFFER_2D_H
