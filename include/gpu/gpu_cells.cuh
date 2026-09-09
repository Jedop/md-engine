#pragma once

__global__
void find_cell_boundaries_kernel(const int* cell_ids, int* cell_start, int* cell_end, int N);

__global__
void calculate_cell_ids_kernel(const double* pos_x, const double* pos_y, const double* pos_z,
                          int* particle_indices, int* cell_ids,
                          int N, double cell_size, int nx);

__global__
void sort_position_kernel(const double* pos_x, const double* pos_y, const double* pos_z,
                        double* sorted_pos_x, double* sorted_pos_y, double* sorted_pos_z,
                        int N, int* particle_indices);