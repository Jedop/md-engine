#include "forces.hpp"

// Computes all forces

std::pair<std::vector<Vec3>, double>
compute_all_forces(const std::vector<Particle> &Particles,
                   const std::vector<int> &head, const std::vector<int> &next,
                   int nx, double cell_size, double box) {

  std::vector<Vec3> all_acc(Particles.size(), {0, 0, 0});

  double potential_energy = 0;
  double inv_rc2 = 1.0 / (rc * rc);
  double inv_rc6 = pow(inv_rc2, 3);
  double U_rc = 4.0 * inv_rc6 * (inv_rc6 - 1.0);
  double box_r = 1 / box;
  // double F_rc = 24 * inv_rc2 * inv_rc6 * (2 * inv_rc6 - 1);

  int num_threads = omp_get_max_threads();

  // Give each thread it's own arrays to avoid race conditions
  std::vector<std::vector<Vec3>> thread_acc(
      num_threads, std::vector<Vec3>(Particles.size(), {0, 0, 0}));
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
#pragma omp for reduction(+ : potential_energy)
    for (size_t i = 0; i < Particles.size(); i++) {

      auto [cx, cy, cz] = get_cell(Particles[i].position, cell_size);

      for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
          for (int dz = -1; dz <= 1; dz++) {
            int nc = cell_index(cx + dx, cy + dy, cz + dz, nx);
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
              thread_acc[tid][i] += F;
              thread_acc[tid][j] += F * -1;
              j = next[j];
            }
          }
        }
      }
    }
  }
#pragma omp parallel for
  for (size_t k = 0; k < Particles.size(); k++) {
    // Each thread gets a chunk of particles (e.g., Thread 0 handles k=0 to
    // 1999) It sums up the forces from ALL threads for its specific particles.
    for (int t = 0; t < num_threads; t++) {
      all_acc[k] += thread_acc[t][k];
    }
  }
  return {all_acc, potential_energy};
}
