#pragma once

#include <cmath>
#include <tuple>
#include <vector>

#include "Particle.hpp"
#include "Vectors.hpp"
#include "constants.hpp"

int cell_index(int ix, int iy, int iz, int nx);

std::tuple<int, int, int> get_cell(const Vec3 &p, double cell_size);

void build_cell_lists(const std::vector<Particle> &Particles,
                      std::vector<int> &head, std::vector<int> &next, int nx,
                      double cell_size);
