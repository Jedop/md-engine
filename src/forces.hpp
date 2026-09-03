#pragma once

#include <cmath>
#include <omp.h>
#include <vector>

#include "Particle.hpp"
#include "Vectors.hpp"
#include "cell_list.hpp"
#include "constants.hpp"

struct GpuMemory {
    double *d_pos_x, *d_pos_y, *d_pos_z;
    double *d_acc_x, *d_acc_y, *d_acc_z;
    double *d_potential_energy;
};

GpuMemory allocate_gpu_memory(int N);
void free_gpu_memory(GpuMemory &mem);

std::pair<std::vector<Vec3>, double>
compute_all_forces(const std::vector<Particle> &Particles,
                   const std::vector<int> &head, const std::vector<int> &next,
                   int nx, double cell_size, double box, GpuMemory &mem);
