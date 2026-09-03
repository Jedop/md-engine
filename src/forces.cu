#include "forces.hpp"
#include <cuda_runtime.h>

__global__
void compute_forces_kernel(const double* pos_x, const double* pos_y, const double* pos_z,
                          double* acc_x, double* acc_y, double* acc_z, double* d_potential_energy,
                          int N, double box, double box_r) {
  
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= N) return;

  const double rc2 = 6.25;              // 2.5^2
  const double inv_rc2 = 1.0 / rc2;
  const double inv_rc6 = inv_rc2 * inv_rc2 * inv_rc2;
  double U_rc = 4.0 * inv_rc6 * (inv_rc6 - 1.0);
  
  double my_x = pos_x[i];
  double my_y = pos_y[i];
  double my_z = pos_z[i];

  double fx = 0.0;
  double fy = 0.0;
  double fz = 0.0;

  double U = 0.0;

  for (int j = 0; j < N; j++) {

    if (i == j) continue;

    double dx = my_x - pos_x[j];
    double dy = my_y - pos_y[j];
    double dz = my_z - pos_z[j];

    dx -= box * round(dx * box_r);
    dy -= box * round(dy * box_r);
    dz -= box * round(dz * box_r);

    double r2 = dx*dx + dy*dy + dz*dz;

    if (r2 < 1e-12 || r2 > 6.25) continue; // 2.5^2 = 6.25

    double inv_r2 = 1.0 / r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;

    // Lennard-Jones Force magnitude
    double f_mag = 24.0 * inv_r2 * inv_r6 * (2.0 * inv_r6 - 1.0);
    U += 0.5 * (4.0 * inv_r6 * (inv_r6 - 1.0) - U_rc);

    fx += dx * f_mag;
    fy += dy * f_mag;
    fz += dz * f_mag;
  }

  acc_x[i] = fx;
  acc_y[i] = fy;
  acc_z[i] = fz;
  d_potential_energy[i] = U;
}
// Computes all forces

std::pair<std::vector<Vec3>, double>
compute_all_forces(const std::vector<Particle> &Particles,
                   const std::vector<int> &head, const std::vector<int> &next,
                   int nx, double cell_size, double box, GpuMemory &mem) {
   int N = Particles.size();
   double box_r = 1 / box; 
   double potential_energy = 0;

   // 1. Convert AoS to SoA (C++ Structs to Flat Arrays)
    std::vector<double> h_pos_x(N), h_pos_y(N), h_pos_z(N);
    for(int i=0; i<N; i++) {
        h_pos_x[i] = Particles[i].position.x;
        h_pos_y[i] = Particles[i].position.y;
        h_pos_z[i] = Particles[i].position.z;
    }

    size_t bytes = N * sizeof(double);
    
    // 2. Copy Data to GPU
    cudaMemcpy(mem.d_pos_x, h_pos_x.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(mem.d_pos_y, h_pos_y.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(mem.d_pos_z, h_pos_z.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemset(mem.d_acc_x, 0, bytes); 
    cudaMemset(mem.d_acc_y, 0, bytes); 
    cudaMemset(mem.d_acc_z, 0, bytes);
    cudaMemset(mem.d_potential_energy, 0, bytes);
    
    // 4. Launch Kernel
    int threads = 256;
    int blocks = (N + threads - 1) / threads; // Ceiling division
    compute_forces_kernel<<<blocks, threads>>>(mem.d_pos_x, mem.d_pos_y, mem.d_pos_z, mem.d_acc_x, mem.d_acc_y, mem.d_acc_z, mem.d_potential_energy, N, box, box_r);
    
    // Wait for GPU to finish (Good for debugging)
    cudaDeviceSynchronize(); 

    // 5. Copy Data back to CPU
    std::vector<double> h_acc_x(N), h_acc_y(N), h_acc_z(N), h_pot(N);
    cudaMemcpy(h_acc_x.data(), mem.d_acc_x, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_acc_y.data(), mem.d_acc_y, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_acc_z.data(), mem.d_acc_z, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_pot.data(), mem.d_potential_energy, bytes, cudaMemcpyDeviceToHost);

    // 7. Repack into C++ Vector and Return
    std::vector<Vec3> all_acc(N);
    for(int i=0; i<N; i++) {
        all_acc[i] = {h_acc_x[i], h_acc_y[i], h_acc_z[i]};
        potential_energy += h_pot[i];
    }

    return {all_acc, potential_energy}; 
   }

GpuMemory allocate_gpu_memory(int N) {
    GpuMemory mem;
    size_t bytes = N * sizeof(double);
    
    cudaMalloc(&mem.d_pos_x, bytes);
    cudaMalloc(&mem.d_pos_y, bytes);
    cudaMalloc(&mem.d_pos_z, bytes);
    cudaMalloc(&mem.d_acc_x, bytes);
    cudaMalloc(&mem.d_acc_y, bytes);
    cudaMalloc(&mem.d_acc_z, bytes);
    cudaMalloc(&mem.d_potential_energy, bytes);
    
    return mem;
}

void free_gpu_memory(GpuMemory &mem) {
    cudaFree(mem.d_pos_x); cudaFree(mem.d_pos_y); cudaFree(mem.d_pos_z);
    cudaFree(mem.d_acc_x); cudaFree(mem.d_acc_y); cudaFree(mem.d_acc_z);
    cudaFree(mem.d_potential_energy);
}

