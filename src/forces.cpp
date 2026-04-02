#include "forces.hpp"

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
