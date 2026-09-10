#pragma once

#include <vector>
#include "Particle.hpp"

struct GpuMemory {
    double *d_pos_x, *d_pos_y, *d_pos_z;
    float *d_sorted_pos_x, *d_sorted_pos_y, *d_sorted_pos_z;
    double *d_vel_x, *d_vel_y, *d_vel_z;
    double *d_acc_x, *d_acc_y, *d_acc_z;
    double *d_potential_energy;
    double *d_kinetic_energy;

    int *d_particle_indices;
    int *d_cell_ids;

    int *d_cell_start;
    int *d_cell_end;
};

GpuMemory allocate_gpu_memory(int N, int num_cells);
void free_gpu_memory(GpuMemory &mem);
void init_gpu_memory(std::vector<Particle> &Particles, GpuMemory &mem);
void sync_gpu_to_cpu(std::vector<Particle> &Particles, GpuMemory &mem);