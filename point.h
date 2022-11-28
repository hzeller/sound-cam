#ifndef POINT_H_
#define POINT_H_

#include <math.h>

class Point {
public:
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
  
private:
  static float square(float a) { return a * a; }
};

constexpr Point operator+(const Point &a, const Point &b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
constexpr Point operator-(const Point &a, const Point &b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
constexpr Point operator*(const Point &p, float scalar) {
  return {p.x * scalar, p.y * scalar, p.z * scalar};
}

#endif // POINT_H_
