#include "cell_list.hpp"

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
