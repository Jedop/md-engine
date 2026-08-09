#include "thermostat.hpp"

void apply_velocity_rescaling(std::vector<Particle> &Particles,
                              double T_current, double T_target, double dt) {
  // Literally just correct the temperature in one step                                
  double lambda = std::sqrt(T_target / T_current);
#pragma omp parallel for
  for (auto &p : Particles) {
    p.velocity = p.velocity * lambda;
  }
}

void apply_berendsen_thermostat(std::vector<Particle> &Particles,
                                double T_current, double T_target, double dt) {
  // tau is the coupling time.
  double tau = 0.1;

  // The Berendsen scaling factor
  double ratio = T_target / T_current;
  double lambda = std::sqrt(1.0 + (dt / tau) * (ratio - 1.0));

// Apply the gentle scale
#pragma omp parallel for
  for (auto &p : Particles) {
    p.velocity = p.velocity * lambda;
  }
}
