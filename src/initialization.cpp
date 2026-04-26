#include "initialization.hpp"

std::vector<Particle> init_fcc_lattice(int unit_cells_per_side,
                                       double box_length) {
  std::vector<Particle> Particles;
  double a = box_length / unit_cells_per_side;
  int N = 4 * (unit_cells_per_side * unit_cells_per_side *
               unit_cells_per_side); // 4000 atoms

  for (int i = 0; i < unit_cells_per_side; i++) {
    for (int j = 0; j < unit_cells_per_side; j++) {
      for (int k = 0; k < unit_cells_per_side; k++) {

        // Base corner of the current unit cell
        double x0 = i * a;
        double y0 = j * a;
        double z0 = k * a;

        // The 4 basis atoms of an FCC unit cell
        // 1. Corner atom
        Particles.emplace_back(x0, y0, z0);
        // 2. Face-center XY
        Particles.emplace_back(x0 + a / 2.0, y0 + a / 2.0, z0);
        // 3. Face-center XZ
        Particles.emplace_back(x0 + a / 2.0, y0, z0 + a / 2.0);
        // 4. Face-center YZ
        Particles.emplace_back(x0, y0 + a / 2.0, z0 + a / 2.0);
      }
    }
  }
  return Particles;
}

std::vector<Particle> init_sc_lattice(int particles_per_side,
                                      double box_length) {
  std::vector<Particle> Particles;
  double spacing = box_length / particles_per_side;
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

  return Particles;
}

void remove_center_of_mass_momentum(std::vector<Particle> &Particles) {

  srand(time(0));

  for (size_t i = 0; i < Particles.size(); i++) {

    Particles[i].velocity = {((double)rand()) / RAND_MAX - 0.5,
                             ((double)rand()) / RAND_MAX - 0.5,
                             ((double)rand()) / RAND_MAX - 0.5};
  }

  Vec3 net_p = {0, 0, 0};
  for (const auto &p : Particles)
    net_p += p.velocity;
  Vec3 com_v = net_p / Particles.size();
  for (auto &p : Particles)
    p.velocity = p.velocity - com_v;
}
