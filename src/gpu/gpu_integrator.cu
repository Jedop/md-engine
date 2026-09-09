#include "gpu_integrator.cuh"
#include <cuda_runtime.h>

// Velocity Verlet Step 1
__global__
void verlet_step1_kernel(double* x, double* y, double* z,
                         double* vx, double* vy, double* vz,
                         const double* ax, const double* ay, const double* az,
                         int N, double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    x[i] += vx[i] * dt + 0.5 * ax[i] * dt * dt;
    y[i] += vy[i] * dt + 0.5 * ay[i] * dt * dt;
    z[i] += vz[i] * dt + 0.5 * az[i] * dt * dt;

    vx[i] += 0.5 * ax[i] * dt;
    vy[i] += 0.5 * ay[i] * dt;
    vz[i] += 0.5 * az[i] * dt;
}

// Velocity Verlet Step 2
__global__
void verlet_step2_kernel(double* vx, double* vy, double* vz,
                         const double* ax, const double* ay, const double* az,
                         double* d_kinetic_energy, // Track KE here!
                         int N, double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    vx[i] += 0.5 * ax[i] * dt;
    vy[i] += 0.5 * ay[i] * dt;
    vz[i] += 0.5 * az[i] * dt;

    // Calculate Kinetic Energy
    double v2 = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    d_kinetic_energy[i] = 0.5 * v2;
}