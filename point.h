#ifndef POINT_H_
#define POINT_H_

#include <math.h>

class Point {
public:
  real_t distance_to(const Point &other) const {
    return sqrtf(square(x - other.x) + square(y - other.y) +
                 square(z - other.z));
  }

  real_t length() const { return sqrtf(square(x) + square(y) + square(z)); }

  real_t dotMul(const Point &other) const {
    return x * other.x + y * other.y + z * other.z;
  }

  void MakeUnitLen() {
    real_t len = length();
    x /= len;
    y /= len;
    z /= len;
  }

  real_t x, y, z;

private:
  static real_t square(real_t a) { return a * a; }
};

constexpr Point operator+(const Point &a, const Point &b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
constexpr Point operator-(const Point &a, const Point &b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
constexpr Point operator*(const Point &p, real_t scalar) {
  return {p.x * scalar, p.y * scalar, p.z * scalar};
}

#endif // POINT_H_
