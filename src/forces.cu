#include "forces.hpp"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

__global__
void find_cell_boundaries_kernel(const int* cell_ids, int* cell_start, int* cell_end, int N) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    int my_cell = cell_ids[i];

    if (i == 0) {
        cell_start[my_cell] = i;
    } else if (my_cell != cell_ids[i - 1]) {
        cell_start[my_cell] = i;
    }

    if (i == N - 1) {
        cell_end[my_cell] = N;
    } else if (my_cell != cell_ids[i + 1]) {
        cell_end[my_cell] = i + 1;
    }
}

__global__
void calculate_cell_ids_kernel(const double* pos_x, const double* pos_y, const double* pos_z,
                          int* particle_indices, int* cell_ids,
                          int N, double cell_size, int nx) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N) return;

    int ix = (int(floor(pos_x[i] / cell_size)) % nx + nx) % nx;
    int iy = (int(floor(pos_y[i] / cell_size)) % nx + nx) % nx;
    int iz = (int(floor(pos_z[i] / cell_size)) % nx + nx) % nx;

    int c_id = ix + iy * nx + iz * nx * nx;

    particle_indices[i] = i;
    cell_ids[i] = c_id;
}

__global__
void sort_position_kernel(const double* pos_x, const double* pos_y, const double* pos_z,
                        double* sorted_pos_x, double* sorted_pos_y, double* sorted_pos_z,
                        int N, int* particle_indices) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= N) return;
    
    sorted_pos_x[i] = pos_x[particle_indices[i]];
    sorted_pos_y[i] = pos_y[particle_indices[i]];
    sorted_pos_z[i] = pos_z[particle_indices[i]];
}

__global__
void compute_forces_kernel(const double* pos_x, const double* pos_y, const double* pos_z,
                          double* acc_x, double* acc_y, double* acc_z, double* d_potential_energy,
                          int* particle_indices, int* cell_ids, int* cell_start, int* cell_end,
                          int N, double box, double box_r, int nx) {
  
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N) return;

    const double rc2 = 6.25;              // 2.5^2
    const double inv_rc2 = 1.0 / rc2;
    const double inv_rc6 = inv_rc2 * inv_rc2 * inv_rc2;
    double U_rc = 4.0 * inv_rc6 * (inv_rc6 - 1.0);

    double my_x = pos_x[i];
    double my_y = pos_y[i];
    double my_z = pos_z[i];
    int my_cell = cell_ids[i];

    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;

    double U = 0.0;

    int cz = my_cell / (nx * nx);
    int cy = (my_cell / nx) % nx;
    int cx = my_cell % nx;
    
    for (int dz = -1; dz <= 1; dz++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                int n_x = (cx + dx + nx) % nx;
                int n_y = (cy + dy + nx) % nx;
                int n_z = (cz + dz + nx) % nx;

                int neighbor_cell = n_x + n_y * nx + n_z * nx * nx;

                int start = cell_start[neighbor_cell];
                int end = cell_end[neighbor_cell];

                if (start == -1) continue;

                for (int j = start; j < end; j++) {
                    if (i == j) continue; // Don't interact with yourself
                    
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
            }
        }
    }


    int original_idx = particle_indices[i];


    acc_x[original_idx] = fx;
    acc_y[original_idx] = fy;
    acc_z[original_idx] = fz;
    d_potential_energy[original_idx] = U;
    }
// Computes all forces

std::pair<std::vector<Vec3>, double>
compute_all_forces_gpu(const std::vector<Particle> &Particles,
                        double box, GpuMemory &mem) {
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
    const int nx = int(box / rc);
    const double cell_size = box / nx;
    const int num_cells = nx * nx * nx;

    calculate_cell_ids_kernel<<<blocks, threads>>>(mem.d_pos_x, mem.d_pos_y, mem.d_pos_z, mem.d_particle_indices, mem.d_cell_ids, N, cell_size, nx);

    cudaDeviceSynchronize();
    
    thrust::sort_by_key(thrust::device, 
                        mem.d_cell_ids, mem.d_cell_ids + N, 
                        mem.d_particle_indices);

    sort_position_kernel<<<blocks, threads>>>(mem.d_pos_x, mem.d_pos_y, mem.d_pos_z, mem.d_sorted_pos_x, mem.d_sorted_pos_y, mem.d_sorted_pos_z, N, mem.d_particle_indices);
    
    cudaDeviceSynchronize();
    
    cudaMemset(mem.d_cell_start, 0xff, num_cells * sizeof(int));
    cudaMemset(mem.d_cell_end, 0xff, num_cells * sizeof(int));

    find_cell_boundaries_kernel<<<blocks, threads>>>(mem.d_cell_ids, 
                                                     mem.d_cell_start, mem.d_cell_end, 
                                                     N);

    cudaDeviceSynchronize();

    compute_forces_kernel<<<blocks, threads>>>(mem.d_sorted_pos_x, mem.d_sorted_pos_y, mem.d_sorted_pos_z, 
        mem.d_acc_x, mem.d_acc_y, mem.d_acc_z, 
        mem.d_potential_energy, mem.d_particle_indices, mem.d_cell_ids, mem.d_cell_start, mem.d_cell_end,
        N, box, box_r, nx);
    
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

GpuMemory allocate_gpu_memory(int N, int num_cells) {
    GpuMemory mem;
    size_t bytes = N * sizeof(double);
    size_t int_bytes = N * sizeof(int);
    size_t cell_bytes = num_cells * sizeof(int);

    cudaMalloc(&mem.d_pos_x, bytes);
    cudaMalloc(&mem.d_pos_y, bytes);
    cudaMalloc(&mem.d_pos_z, bytes);
    cudaMalloc(&mem.d_sorted_pos_x, bytes);
    cudaMalloc(&mem.d_sorted_pos_y, bytes);
    cudaMalloc(&mem.d_sorted_pos_z, bytes);
    cudaMalloc(&mem.d_acc_x, bytes);
    cudaMalloc(&mem.d_acc_y, bytes);
    cudaMalloc(&mem.d_acc_z, bytes);
    cudaMalloc(&mem.d_potential_energy, bytes);
    cudaMalloc(&mem.d_particle_indices, int_bytes);
    cudaMalloc(&mem.d_cell_ids, int_bytes);
    cudaMalloc(&mem.d_cell_start, cell_bytes);
    cudaMalloc(&mem.d_cell_end, cell_bytes);
    
    return mem;
}

void free_gpu_memory(GpuMemory &mem) {
    cudaFree(mem.d_pos_x); cudaFree(mem.d_pos_y); cudaFree(mem.d_pos_z);
    cudaFree(mem.d_sorted_pos_x); cudaFree(mem.d_sorted_pos_y); cudaFree(mem.d_sorted_pos_z);
    cudaFree(mem.d_acc_x); cudaFree(mem.d_acc_y); cudaFree(mem.d_acc_z);
    cudaFree(mem.d_potential_energy);

    cudaFree(mem.d_particle_indices);
    cudaFree(mem.d_cell_ids);
    cudaFree(mem.d_cell_start);
    cudaFree(mem.d_cell_end);
}

