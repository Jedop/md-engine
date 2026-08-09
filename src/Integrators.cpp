#include "Integrators.hpp"

// Euler
void Euler::step1(std::vector<Particle> &Particles, double dt) {
  for (auto &p : Particles) {
    p.position += p.velocity * dt;
    p.velocity += p.acceleration * dt;
  }
}

void Euler::step2(std::vector<Particle> &Particles,
                  const std::vector<Vec3> &new_accelerations, double dt) {
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].acceleration = new_accelerations[i];
  }
}

// Leapfrog
void Leapfrog::step1(std::vector<Particle> &Particles, double dt) {
  for (auto &p : Particles) {
    p.position += p.velocity * (dt / 2.0);
  }
}

void Leapfrog::step2(std::vector<Particle> &Particles,
                     const std::vector<Vec3> &new_accelerations, double dt) {
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].acceleration = new_accelerations[i];
    Particles[i].velocity += Particles[i].acceleration * dt;
    Particles[i].position += Particles[i].velocity * (dt / 2.0);
  }
}

// VelocityVerlet
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
