#pragma once
#include "Vectors.hpp"

struct Particle {
  Vec3 position;
  Vec3 velocity;
  Vec3 acceleration;

  Particle(double x, double y, double z) {
    position = {x, y, z};
    velocity = {0, 0, 0};
    acceleration = {0, 0, 0};
  }
};
