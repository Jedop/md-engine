#include "simulation.hpp"

std::tuple<double, double, double> update(std::vector<Particle> &Particles,
                                          std::vector<int> &head,
                                          std::vector<int> &next, double dt) {
  // 1. Pre-Force Update (Move everyone to new positions)
  VelocityVerlet::step1(Particles, dt);

  // 2. Rebuild the Grid
  build_cell_lists(Particles, head, next);

  // 3. Synchronous Force Calculation (Calculate forces at the NEW positions)
  auto [new_acc, potential_energy] = compute_all_forces(Particles, head, next);

  // 4. Post-Force Update (Finish the time step using the new forces)
  VelocityVerlet::step2(Particles, new_acc, dt);

  // Calculates Kinetic Energy
  double kinetic_energy = 0;

  for (const auto &p : Particles) {
    Vec3 v = p.velocity;
    // v = p.velocity - p.acceleration * dt * 0.5; Leapfrog Calculation
    double v2 = v.x * v.x + v.y * v.y + v.z * v.z;
    kinetic_energy += 0.5 * v2;
  }

  double Energy = potential_energy + kinetic_energy;
  return {potential_energy, kinetic_energy, Energy};
}
