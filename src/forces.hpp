#pragma once

#include <cmath>
#include <omp.h>
#include <vector>

#include "Particle.hpp"
#include "Vectors.hpp"
#include "cell_list.hpp"
#include "constants.hpp"

std::pair<std::vector<Vec3>, double>
compute_all_forces(const std::vector<Particle> &Particles,
                   const std::vector<int> &head, const std::vector<int> &next,
                   int nx, double cell_size, double box);
