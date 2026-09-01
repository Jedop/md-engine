#include "forces.hpp"
#include <cuda_runtime.h>
#include <iostream>

__global__
void compute_forces_kernel(const double* pos_x, const double* pos_y, const double* pos_z,
                          double* acc_x, double* acc_y, double* acc_z,
                          int N, double box, double box_r) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= N) return;
  
  float my_x = pos_x[i];
  float my_y = pos_y[i];
  double my_z = pos_z[i];

  double fx = 0.0f;
  double fy = 0.0f;
  double fz = 0.0f;

  for (int j = 0; j < N; j++) {

    if (i == j) continue;

    double dx = my_x - pos_x[j];
    double dy = my_y - pos_y[j];
    double dz = my_z - pos_z[j];

    dx -= box * roundf(dx * box_r);
    dy -= box * roundf(dy * box_r);
    dz -= box * roundf(dz * box_r);

    double r2 = dx*dx + dy*dy + dz*dz;

    if (r2 < 1e-12f || r2 > 6.25f) continue; // 2.5^2 = 6.25

    double inv_r2 = 1.0f / r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;

    // Lennard-Jones Force magnitude
    double f_mag = 24.0f * inv_r2 * inv_r6 * (2.0f * inv_r6 - 1.0f);

    fx += dx * f_mag;
    fy += dy * f_mag;
    fz += dz * f_mag;
  }

  acc_x[i] = fx;
  acc_y[i] = fy;
  acc_z[i] = fz;
}
// Computes all forces

std::pair<std::vector<Vec3>, double>
compute_all_forces(const std::vector<Particle> &Particles,
                   const std::vector<int> &head, const std::vector<int> &next,
                   int nx, double cell_size, double box) {
   int N = Particles.size();
   double box_r = 1 / box; 

   // 1. Convert AoS to SoA (C++ Structs to Flat Arrays)
    std::vector<double> h_pos_x(N), h_pos_y(N), h_pos_z(N);
    for(int i=0; i<N; i++) {
        h_pos_x[i] = Particles[i].position.x;
        h_pos_y[i] = Particles[i].position.y;
        h_pos_z[i] = Particles[i].position.z;
    }

    // 2. Allocate GPU Memory (Device)
    double *d_pos_x, *d_pos_y, *d_pos_z;
    double *d_acc_x, *d_acc_y, *d_acc_z;
    size_t bytes = N * sizeof(double);
    
    cudaMalloc(&d_pos_x, bytes); cudaMalloc(&d_pos_y, bytes); cudaMalloc(&d_pos_z, bytes);
    cudaMalloc(&d_acc_x, bytes); cudaMalloc(&d_acc_y, bytes); cudaMalloc(&d_acc_z, bytes);

    // 3. Copy Data to GPU
    cudaMemcpy(d_pos_x, h_pos_x.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_pos_y, h_pos_y.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_pos_z, h_pos_z.data(), bytes, cudaMemcpyHostToDevice);
    // Zero out the acceleration arrays on the GPU
    cudaMemset(d_acc_x, 0, bytes); cudaMemset(d_acc_y, 0, bytes); cudaMemset(d_acc_z, 0, bytes);

    // 4. Launch Kernel
    int threads = 256;
    int blocks = (N + threads - 1) / threads; // Ceiling division
    compute_forces_kernel<<<blocks, threads>>>(d_pos_x, d_pos_y, d_pos_z, d_acc_x, d_acc_y, d_acc_z, N, box, box_r);
    
    // Wait for GPU to finish (Good for debugging)
    cudaDeviceSynchronize(); 

    // 5. Copy Data back to CPU
    std::vector<double> h_acc_x(N), h_acc_y(N), h_acc_z(N);
    cudaMemcpy(h_acc_x.data(), d_acc_x, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_acc_y.data(), d_acc_y, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_acc_z.data(), d_acc_z, bytes, cudaMemcpyDeviceToHost);

    // 6. Free GPU Memory
    cudaFree(d_pos_x); cudaFree(d_pos_y); cudaFree(d_pos_z);
    cudaFree(d_acc_x); cudaFree(d_acc_y); cudaFree(d_acc_z);

    // 7. Repack into C++ Vector and Return
    std::vector<Vec3> all_acc(N);
    for(int i=0; i<N; i++) {
        all_acc[i] = {h_acc_x[i], h_acc_y[i], h_acc_z[i]};
    }

    return {all_acc, 0.0}; 
   }

