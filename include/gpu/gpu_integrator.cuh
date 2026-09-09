#pragma once

#include "gpu_memory.hpp"

__global__
void verlet_step1_kernel(double* x, double* y, double* z,
                         double* vx, double* vy, double* vz,
                         const double* ax, const double* ay, const double* az,
                         int N, double dt);

__global__
void verlet_step2_kernel(double* vx, double* vy, double* vz,
                         const double* ax, const double* ay, const double* az,
                         double* d_kinetic_energy, // Track KE here!
                         int N, double dt);

