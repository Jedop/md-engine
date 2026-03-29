#include "Integrators.hpp"

void VelocityVerlet::step1(std::vector<Particle> &Particles, double dt) {
  for (auto &p : Particles) {
    p.position += p.velocity * dt + p.acceleration * (dt * dt / 2.0);
    p.velocity += p.acceleration * (dt / 2.0);
  }
}

void VelocityVerlet::step2(std::vector<Particle> &Particles,
                           const std::vector<Vec3> &new_accelerations,
                           double dt) {
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].acceleration = new_accelerations[i];
    Particles[i].velocity += Particles[i].acceleration * (dt / 2.0);
  }
}
