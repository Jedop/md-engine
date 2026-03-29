#pragma once

struct Vec3 {
  double x, y, z;

  Vec3 operator+(const Vec3 &o) const { return {x + o.x, y + o.y, z + o.z}; }

  Vec3 operator-(const Vec3 &o) const { return {x - o.x, y - o.y, z - o.z}; }

  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }

  Vec3 &operator+=(const Vec3 &o) {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }

  void apply_pbc(double box, double box_r) {
    auto wrap = [&](double &d) {
      int k = (int)(d * box_r + (d >= 0 ? 0.5 : -0.5));
      d -= k * box;
    };
    wrap(x);
    wrap(y);
    wrap(z);
  }
};
