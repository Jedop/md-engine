#include "gpu_cells.cuh"
#include <cuda_runtime.h>

// Finds the first and last elements of each cell, and stores them
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

// Calculates which cell each particle belongs to
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

// Sorts the positions according to the particle indices, so we can do spatial hashing to calculate the forces easily
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