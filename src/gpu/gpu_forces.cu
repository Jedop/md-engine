#include "gpu_forces.hpp"
#include "gpu_cells.cuh"
#include "gpu_integrator.cuh"
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

// Calculates the forces on the GPU
__global__
void compute_forces_kernel(const double* pos_x, const double* pos_y, const double* pos_z,
                          double* acc_x, double* acc_y, double* acc_z, double* d_potential_energy,
                          int* particle_indices, int* cell_ids, int* cell_start, int* cell_end,
                          int N, double box, double box_r, int nx) {
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N) return;

    // Ideally should somehow use the constant from constants.hpp, but I am too lazy to figure that out at this point. 
    // Probably easiest to simply pass it as a parameter to the function
    const float rc2 = 6.25;              // 2.5^2
    const float inv_rc2 = 1.0 / rc2;
    const float inv_rc6 = inv_rc2 * inv_rc2 * inv_rc2;
    const float U_rc = 4.0 * inv_rc6 * (inv_rc6 - 1.0);

    const float box_f = (float)box;
    const float box_r_f = (float)box_r;
    
    // Position of ith particle
    float my_x = (float)pos_x[i];
    float my_y = (float)pos_y[i];
    float my_z = (float)pos_z[i];
    int my_cell = cell_ids[i];
    
    // Magnitude of force on ith particle
    float fx = 0.0f;
    float fy = 0.0f;
    float fz = 0.0f;
    float U = 0.0f;
    
    // Cell index
    int cz = my_cell / (nx * nx);
    int cy = (my_cell / nx) % nx;
    int cx = my_cell % nx;
    
    for (int dz = -1; dz <= 1; dz++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                
                // Adjacent cells
                int n_x = (cx + dx + nx) % nx;
                int n_y = (cy + dy + nx) % nx;
                int n_z = (cz + dz + nx) % nx;

                int neighbor_cell = n_x + n_y * nx + n_z * nx * nx;

                // Find all the particles in that cell
                int start = cell_start[neighbor_cell];
                int end = cell_end[neighbor_cell];

                // start == -1 if there are no particles in that cell
                if (start == -1) continue;

                for (int j = start; j < end; j++) {
                    if (i == j) continue; // Don't interact with yourself

                    // Difference in positions of particles
                    float ddx = my_x - (float)pos_x[j];
                    float ddy = my_y - (float)pos_y[j];
                    float ddz = my_z - (float)pos_z[j];

                    // Minimum Image Convention
                    ddx -= box_f * roundf(ddx * box_r_f);
                    ddy -= box_f * roundf(ddy * box_r_f);
                    ddz -= box_f * roundf(ddz * box_r_f);

                    // r^2 where r is vector pointing from particle 1 to 2
                    float r2 = ddx*ddx + ddy*ddy + ddz*ddz;

                    // Minimum Cutoff and Maximum Cutoff
                    if (r2 < 1e-6f || r2 > 6.25f) continue;

                    // Following is just arithmetic
                    float inv_r2 = 1.0f / r2;
                    float inv_r6 = inv_r2 * inv_r2 * inv_r2;

                    float f_mag = 24.0f * inv_r2 * inv_r6 * (2.0f * inv_r6 - 1.0f);
                    U += 0.5f * (4.0f * inv_r6 * (inv_r6 - 1.0f) - U_rc);

                    fx += ddx * f_mag;
                    fy += ddy * f_mag;
                    fz += ddz * f_mag;
                }
            }
        }
    }
    // Assign the accelerations
    int original_idx = particle_indices[i];

    acc_x[original_idx] = (double)fx;
    acc_y[original_idx] = (double)fy;
    acc_z[original_idx] = (double)fz;
    d_potential_energy[i] = (double)U;
    }

// Computes all forces without the velocity steps, for the stuff before the first step
void init_gpu_physics(GpuMemory &mem, int N, double box, double box_r, int nx, double cell_size) {
    int threads = 256;
    int blocks = (N + threads - 1) / threads;
    int num_cells = nx * nx * nx;

    // Hash, Sort, Boundary
    calculate_cell_ids_kernel<<<blocks, threads>>>(mem.d_pos_x, mem.d_pos_y, mem.d_pos_z, mem.d_particle_indices, mem.d_cell_ids, N, cell_size, nx);
    thrust::sort_by_key(thrust::device, mem.d_cell_ids, mem.d_cell_ids + N, mem.d_particle_indices);
    sort_position_kernel<<<blocks, threads>>>(mem.d_pos_x, mem.d_pos_y, mem.d_pos_z, mem.d_sorted_pos_x, mem.d_sorted_pos_y, mem.d_sorted_pos_z, N, mem.d_particle_indices);
    cudaMemset(mem.d_cell_start, 0xff, num_cells * sizeof(int));
    cudaMemset(mem.d_cell_end, 0xff, num_cells * sizeof(int));
    find_cell_boundaries_kernel<<<blocks, threads>>>(mem.d_cell_ids, mem.d_cell_start, mem.d_cell_end, N);
    
    // Compute initial forces
    cudaMemset(mem.d_acc_x, 0, N * sizeof(double)); 
    cudaMemset(mem.d_acc_y, 0, N * sizeof(double)); 
    cudaMemset(mem.d_acc_z, 0, N * sizeof(double));
    cudaMemset(mem.d_potential_energy, 0, N * sizeof(double));

    compute_forces_kernel<<<blocks, threads>>>(mem.d_sorted_pos_x, mem.d_sorted_pos_y, mem.d_sorted_pos_z, 
        mem.d_acc_x, mem.d_acc_y, mem.d_acc_z, mem.d_potential_energy, 
        mem.d_particle_indices, mem.d_cell_ids, mem.d_cell_start, mem.d_cell_end, N, box, box_r, nx);
}

std::pair<double, double>
compute_all_forces_gpu(const std::vector<Particle> &Particles,
                        double box, GpuMemory &mem, double dt) {
   int N = Particles.size();
   double box_r = 1 / box; 

    int threads = 256;
    int blocks = (N + threads - 1) / threads; // Ceiling division
    const int nx = int(box / rc);
    const double cell_size = box / nx;
    const int num_cells = nx * nx * nx;
    
    // Set Potential and Kinetic Energy to 0 since they need to be reset since they are only added to in the functions
    cudaMemset(mem.d_potential_energy, 0, N * sizeof(double));
    cudaMemset(mem.d_kinetic_energy, 0, N *sizeof(double));
    
    // Verlet Velocity 1st step but as a cuda kernel
    verlet_step1_kernel<<<blocks, threads>>>(mem.d_pos_x, mem.d_pos_y, mem.d_pos_z,
                                        mem.d_vel_x, mem.d_vel_y, mem.d_vel_z,
                                        mem.d_acc_x, mem.d_acc_y, mem.d_acc_z,
                                        N, dt);                            
    
    // Calculate the cell ids
    calculate_cell_ids_kernel<<<blocks, threads>>>(mem.d_pos_x, mem.d_pos_y, mem.d_pos_z, mem.d_particle_indices, mem.d_cell_ids, N, cell_size, nx);


    // Sort the particles wrt to their cell ids
    thrust::sort_by_key(thrust::device, 
                        mem.d_cell_ids, mem.d_cell_ids + N, 
                        mem.d_particle_indices);
    
    // Sort the sorted position arrays
    sort_position_kernel<<<blocks, threads>>>(mem.d_pos_x, mem.d_pos_y, mem.d_pos_z, mem.d_sorted_pos_x, mem.d_sorted_pos_y, mem.d_sorted_pos_z, N, mem.d_particle_indices);

    // These need to be reset outside the kernels as well(could be shifted inside)
    cudaMemset(mem.d_cell_start, 0xff, num_cells * sizeof(int));
    cudaMemset(mem.d_cell_end, 0xff, num_cells * sizeof(int));
    
    // Self-explanatory
    find_cell_boundaries_kernel<<<blocks, threads>>>(mem.d_cell_ids, 
                                                     mem.d_cell_start, mem.d_cell_end, 
                                                     N);
    
    // The main kernel which computes the forces                                                    
    compute_forces_kernel<<<blocks, threads>>>(mem.d_sorted_pos_x, mem.d_sorted_pos_y, mem.d_sorted_pos_z, 
        mem.d_acc_x, mem.d_acc_y, mem.d_acc_z, 
        mem.d_potential_energy, mem.d_particle_indices, mem.d_cell_ids, mem.d_cell_start, mem.d_cell_end,
        N, box, box_r, nx);
    
    // Velocity Verlet step 2 but as a cuda kernel
    verlet_step2_kernel<<<blocks, threads>>>(mem.d_vel_x, mem.d_vel_y, mem.d_vel_z,
                                        mem.d_acc_x, mem.d_acc_y, mem.d_acc_z,
                                        mem.d_kinetic_energy, N, dt);

    // Synchronize
    cudaDeviceSynchronize();
    
    // Convert the array into a double
    // Thrust automatically performs a highly optimized parallel reduction tree
    thrust::device_ptr<double> pot_ptr(mem.d_potential_energy);
    thrust::device_ptr<double> kin_ptr(mem.d_kinetic_energy);

    double potential_energy = thrust::reduce(pot_ptr, pot_ptr + N, 0.0);
    double kinetic_energy = thrust::reduce(kin_ptr, kin_ptr + N, 0.0);

    return {potential_energy, kinetic_energy};
   }
