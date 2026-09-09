#include "gpu_thermostat.hpp"
#include <cuda_runtime.h>

// A kernel which simply multiplies the velocity by lambda
__global__
void apply_thermostat_kernel(int N, double* vel_x, double* vel_y, double* vel_z, double lambda) {
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < N) {
        vel_x[i] *= lambda;
        vel_y[i] *= lambda;
        vel_z[i] *= lambda;
    }
}

// A kernel which reverses the velocities, thought this was the most appropriate and convenient place to put it
// Arguably gpu_integrator is better, but that causes with header files and I thought that was not worth fixing
__global__
void reverse_velocities_kernel(int N, double* vel_x, double* vel_y, double* vel_z) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < N) {
        vel_x[i] *= -1.0;
        vel_y[i] *= -1.0;
        vel_z[i] *= -1.0;
    }
}

// Function which does what a berendsen thermostat should
void apply_berendsen_thermostat_gpu(std::vector<Particle> &Particles,
                                double T_current, double T_target, double dt, GpuMemory &mem) {
  // tau is the coupling time.
  double tau = 0.1;
  int N = Particles.size();

  // The Berendsen scaling factor
  double ratio = T_target / T_current;
  double lambda = std::sqrt(1.0 + (dt / tau) * (ratio - 1.0));

//   std::cout << "T_current = " << T_current
//           << ", T_target = " << T_target
//           << ", lambda = " << lambda << "\n";

  int threads = 256;
  int blocks = (N + threads - 1) / threads;

  apply_thermostat_kernel<<<blocks, threads>>>(N, mem.d_vel_x, mem.d_vel_y, mem.d_vel_z, lambda);
}

// Reverses velocities, simply a wrapper to call the kernel
void reverse_velocities_gpu(int N, GpuMemory& mem) {
    constexpr int threads = 256;
    int blocks = (N + threads - 1) / threads;

    reverse_velocities_kernel<<<blocks, threads>>>(
        N,
        mem.d_vel_x,
        mem.d_vel_y,
        mem.d_vel_z
    );

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        throw std::runtime_error(cudaGetErrorString(err));
    }
}