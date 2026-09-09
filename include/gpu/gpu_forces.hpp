#pragma once

#include <cmath>
#include <omp.h>
#include <vector>

#include "Particle.hpp"
#include "Vectors.hpp"
#include "cell_list.hpp"
#include "constants.hpp"
#include "gpu_memory.hpp"

std::pair<double, double>
compute_all_forces_gpu(const std::vector<Particle> &Particles,
                    double box, GpuMemory &mem, double dt);

std::pair<std::vector<Vec3>, double>
compute_all_forces(const std::vector<Particle> &Particles,
                   const std::vector<int> &head, const std::vector<int> &next,
                   int nx, double cell_size, double box);

void init_gpu_physics(GpuMemory &mem, int N, double box, double box_r, int nx, double cell_size);