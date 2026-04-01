#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include "Integrators.hpp"
#include "Particle.hpp"
#include "Vectors.hpp"

constexpr int no_of_frames = 10000;
constexpr double dt = 0.001;
constexpr int particles_per_side = 20;
constexpr int N = particles_per_side * particles_per_side * particles_per_side;
constexpr double box = particles_per_side * 1.122;
constexpr double box_r = 1 / box;
constexpr double rc = 2.5;

constexpr int nx = int(box / rc);
constexpr int num_cells = nx * nx * nx;
constexpr double cell_size = box / nx;

// Calculates flat index of cell
int cell_index(int ix, int iy, int iz) {
  ix = (ix + nx) % nx;
  iy = (iy + nx) % nx;
  iz = (iz + nx) % nx;
  return ix + iy * nx + iz * nx * nx;
}

// Calculates cell index given position
std::tuple<int, int, int> get_cell(const Vec3 &p) {
  int ix = int(std::floor(p.x / cell_size));
  int iy = int(std::floor(p.y / cell_size));
  int iz = int(std::floor(p.z / cell_size));
  return {ix, iy, iz};
}

// Builds the cell linked lists
void build_cell_lists(const std::vector<Particle> &Particles,
                      std::vector<int> &head, std::vector<int> &next) {
  std::fill(head.begin(), head.end(), -1);

  for (size_t i = 0; i < Particles.size(); i++) {
    auto [cx, cy, cz] = get_cell(Particles[i].position);
    int c = cell_index(cx, cy, cz);
    next[i] = head[c];
    head[c] = i;
  }
}

// Computes all forces
std::pair<std::vector<Vec3>, double>
compute_all_forces(const std::vector<Particle> &Particles,
                   const std::vector<int> &head, const std::vector<int> &next) {
  std::vector<Vec3> all_acc(Particles.size(), {0, 0, 0});
  double potential_energy = 0;

  double inv_rc2 = 1.0 / (rc * rc);
  double inv_rc6 = pow(inv_rc2, 3);

  double U_rc = 4.0 * inv_rc6 * (inv_rc6 - 1.0);
  // double F_rc = 24 * inv_rc2 * inv_rc6 * (2 * inv_rc6 - 1);

  for (size_t i = 0; i < Particles.size(); i++) {
    auto [cx, cy, cz] = get_cell(Particles[i].position);

    for (int dx = -1; dx <= 1; dx++) {
      for (int dy = -1; dy <= 1; dy++) {
        for (int dz = -1; dz <= 1; dz++) {
          int nc = cell_index(cx + dx, cy + dy, cz + dz);
          int j = head[nc];

          while (j != -1) {

            if (j <= i) {
              j = next[j];
              continue;
            }
            Vec3 r = Particles[i].position - Particles[j].position;
            r.apply_pbc(box, box_r);
            double r2 = r.x * r.x + r.y * r.y + r.z * r.z; // r^2

            if (r2 < 1e-12 || r2 > rc * rc) {
              j = next[j];
              continue; // Safety
            }
            double inv_r2 = 1.0 / r2;
            double inv_r6 = inv_r2 * inv_r2 * inv_r2;

            Vec3 F = r * (24 * inv_r2 * inv_r6 * (2 * inv_r6 - 1));
            double U = 4.0 * inv_r6 * (inv_r6 - 1.0);

            potential_energy += U - U_rc;
            all_acc[i] += F;
            all_acc[j] += F * -1;
            j = next[j];
          }
        }
      }
    }
  }
  return {all_acc, potential_energy};
}

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
    // Particles[i].velocity += Particles[i].acceleration * dt * 0.5; //
    // Leapfrog Initialization
  }

  // Main Loop
  std::ofstream traj_file("trajectory.xyz");
  std::ofstream data_file("thermo.dat");
  for (int i = 0; i < no_of_frames; i++) {
    std::cout << i << "\n";
    auto [potential_energy, kinetic_energy, Energy] =
        update(Particles, head, next, dt);

    if (i % 100 == 0) {
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
  return 0;
}
