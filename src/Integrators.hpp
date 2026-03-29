#pragma once
#include <vector>

#include "Particle.hpp"
#include "Vectors.hpp"

namespace VelocityVerlet {

void step1(std::vector<Particle> &Particles, double dt);

void step2(std::vector<Particle> &Particles,
           const std::vector<Vec3> &new_accelerations, double dt);

} // namespace VelocityVerlet
