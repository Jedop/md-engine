#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include "Integrators.hpp"
#include "Particle.hpp"
#include "Vectors.hpp"

constexpr double box = 70;
constexpr double box_r = 1 / box;

std::pair<std::vector<Vec3>, double>
compute_all_forces(const std::vector<Particle> &Particles) {
  std::vector<Vec3> all_acc(Particles.size(), {0, 0, 0});
  double potential_energy = 0;

  for (size_t i = 0; i < Particles.size(); i++) {
    for (size_t j = i + 1; j < Particles.size(); j++) {

      Vec3 r = Particles[i].position - Particles[j].position;
      r.apply_pbc(box, box_r);
      double r2 = r.x * r.x + r.y * r.y + r.z * r.z; // r^2
      double rc = 2.5;

      if (r2 < 1e-12 || r2 > rc * rc)
        continue; // Safety

      double inv_r2 = 1.0 / r2;
      double inv_r6 = pow(inv_r2, 3);

      Vec3 F = r * 24 * inv_r2 * inv_r6 * (2 * inv_r6 - 1);

      potential_energy += std::sqrt((F.x * F.x + F.y * F.y + F.z * F.z) * r2);
      all_acc[i] += F;
      all_acc[j] += F * -1;
    }
  }
  return {all_acc, potential_energy};
}

std::tuple<double, double, double> update(std::vector<Particle> &Particles,
                                          double dt) {
  // 1. Pre-Force Update (Move everyone to new positions)
  VelocityVerlet::step1(Particles, dt);

  // 2. Synchronous Force Calculation (Calculate forces at the NEW positions)
  auto [new_acc, potential_energy] = compute_all_forces(Particles);

  // 3. Post-Force Update (Finish the time step using the new forces)
  VelocityVerlet::step2(Particles, new_acc, dt);

  double kinetic_energy = 0;

  for (const auto &p : Particles) {
    double v2 = p.velocity.x * p.velocity.x + p.velocity.y * p.velocity.y +
                p.velocity.z * p.velocity.z;
    kinetic_energy += 0.5 * v2;
  }

  double Energy = potential_energy + kinetic_energy;
  return {potential_energy, kinetic_energy, Energy};
}

int main() {
  std::vector<Particle> Particles;
  Particles.reserve(512);

  int particles_per_side = 8;
  double spacing = box / particles_per_side;

  for (int i = 0; i < particles_per_side; i++) {
    for (int j = 0; j < particles_per_side; j++) {
      for (int k = 0; k < particles_per_side; k++) {

        double x = (i + 0.5) * spacing;
        double y = (j + 0.5) * spacing;
        double z = (k + 0.5) * spacing;

        Particles.emplace_back(x, y, z);
      }
    }
  }

  Particles[3].position = {12.0001, 12.0001, 38.2501};
  Particles[3].velocity = {1.0, 1.0, 1.0};
  auto [initial_acc, _] = compute_all_forces(Particles);
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].acceleration = initial_acc[i];
  }

  int no_of_frames = 100000;
  std::ofstream traj_file("trajectory.xyz");
  std::ofstream data_file("thermo.dat");
  for (int i = 0; i < no_of_frames; i++) {
    std::cout << i << "\n";
    auto [potential_energy, kinetic_energy, Energy] =
        update(Particles, 0.00001);

    if (i % 1 == 0) {
      traj_file << Particles.size() << "\n";

      // THE MAGIC LINE: This tells OVITO you have a cubic box of size L x L x L
      traj_file << "Lattice=\"" << box << " 0.0 0.0 0.0 " << box
                << " 0.0 0.0 0.0 " << box
                << "\" Properties=species:S:1:pos:R:3\n";

      for (const auto &p : Particles) {
        // Output the actual, UNWRAPPED drifting coordinates
        traj_file << "Ar " << p.position.x << " " << p.position.y << " "
                  << p.position.z << "\n";
      }

      data_file << potential_energy << " " << kinetic_energy << " " << Energy
                << " " << (2 * kinetic_energy / 3 * 512) << "\n";
    }
  }
  traj_file.close();
  return 0;
}
