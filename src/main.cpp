#include <cstdlib>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include "Particle.hpp"
#include "Vectors.hpp"
#include "cell_list.hpp"
#include "constants.hpp"
#include "forces.hpp"
#include "simulation.hpp"

int main() {

  // Initialization
  std::vector<Particle> Particles;
  Particles.reserve(N);
  std::vector<int> head(num_cells, -1);
  std::vector<int> next(N);

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

  // Particles[3].position = {12.0001, 12.0001, 38.2501};
  // Particles[3].velocity = {1.0, 1.0, 1.0};
  //
  srand(time(0));

  for (size_t i = 0; i < Particles.size(); i++) {

    Particles[i].velocity = {((double)rand()) / RAND_MAX - 0.5,
                             ((double)rand()) / RAND_MAX - 0.5,
                             ((double)rand()) / RAND_MAX - 0.5};
  }

  build_cell_lists(Particles, head, next);
  auto [initial_acc, _] = compute_all_forces(Particles, head, next);

  for (size_t i = 0; i < Particles.size(); i++) {

    Particles[i].acceleration = initial_acc[i];
    // Leapfrog Initialization
    // Particles[i].velocity += Particles[i].acceleration * dt * 0.5;
  }

  // Main Loop
  std::ofstream traj_file("trajectory.xyz");
  std::ofstream data_file("thermo.dat");

  for (int i = 0; i < no_of_frames; i++) {
    auto [potential_energy, kinetic_energy, Energy] =
        update(Particles, head, next, dt);

    if (i % 100 == 0) {

      std::cout << i << "\n";
      traj_file << Particles.size() << "\n";

      traj_file << "Lattice=\"" << box << " 0.0 0.0 0.0 " << box
                << " 0.0 0.0 0.0 " << box
                << "\" Properties=species:S:1:pos:R:3\n";

      for (const auto &p : Particles) {
        traj_file << "Ar " << p.position.x << " " << p.position.y << " "
                  << p.position.z << "\n";
      }

      data_file << potential_energy << " " << kinetic_energy << " " << Energy
                << " " << (2 * kinetic_energy / (3 * N)) << "\n";
    }
  }
  traj_file.close();
  data_file.close();
  return 0;
}
