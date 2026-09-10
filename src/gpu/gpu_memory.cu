#include "gpu_memory.hpp"
#include <cuda_runtime.h>

// Allocates all the GPU Memory in the GPU, with the pointers initialized in mem
GpuMemory allocate_gpu_memory(int N, int num_cells) {
    GpuMemory mem;
    size_t bytes = N * sizeof(double);
    size_t int_bytes = N * sizeof(int);
    size_t cell_bytes = num_cells * sizeof(int);
    size_t float_bytes = N * sizeof(float);

    cudaMalloc(&mem.d_pos_x, bytes);
    cudaMalloc(&mem.d_pos_y, bytes);
    cudaMalloc(&mem.d_pos_z, bytes);
    cudaMalloc(&mem.d_sorted_pos_x, float_bytes);
    cudaMalloc(&mem.d_sorted_pos_y, float_bytes);
    cudaMalloc(&mem.d_sorted_pos_z, float_bytes);
    cudaMalloc(&mem.d_vel_x, bytes);
    cudaMalloc(&mem.d_vel_y, bytes);
    cudaMalloc(&mem.d_vel_z, bytes);
    cudaMalloc(&mem.d_acc_x, bytes);
    cudaMalloc(&mem.d_acc_y, bytes);
    cudaMalloc(&mem.d_acc_z, bytes);
    cudaMalloc(&mem.d_potential_energy, bytes);
    cudaMalloc(&mem.d_kinetic_energy, bytes);
    cudaMalloc(&mem.d_particle_indices, int_bytes);
    cudaMalloc(&mem.d_cell_ids, int_bytes);
    cudaMalloc(&mem.d_cell_start, cell_bytes);
    cudaMalloc(&mem.d_cell_end, cell_bytes);
    
    return mem;
}

// Sync the GPU memory back to the CPU, used every 100 time steps to write data to files
void sync_gpu_to_cpu(std::vector<Particle> &Particles, GpuMemory &mem) {
    int N = Particles.size();
    size_t bytes = N * sizeof(double);
    std::vector<double> h_x(N), h_y(N), h_z(N), h_vx(N), h_vy(N), h_vz(N), h_ax(N), h_ay(N), h_az(N);

    cudaMemcpy(h_x.data(), mem.d_pos_x, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_y.data(), mem.d_pos_y, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_z.data(), mem.d_pos_z, bytes, cudaMemcpyDeviceToHost);

    cudaMemcpy(h_vx.data(), mem.d_vel_x, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vy.data(), mem.d_vel_y, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vz.data(), mem.d_vel_z, bytes, cudaMemcpyDeviceToHost);

    cudaMemcpy(h_ax.data(), mem.d_acc_x, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_ay.data(), mem.d_acc_y, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_az.data(), mem.d_acc_z, bytes, cudaMemcpyDeviceToHost);

    for(int i=0; i<N; i++) {
        Particles[i].velocity = {h_vx[i], h_vy[i], h_vz[i]};
        Particles[i].position = {h_x[i], h_y[i], h_z[i]};
        Particles[i].acceleration = {h_ax[i], h_ay[i], h_az[i]};
    }
}

// Initialize the GPU memory
void init_gpu_memory(std::vector<Particle> &Particles, GpuMemory &mem) {

    int N = Particles.size();
    // 1. Convert AoS to SoA (C++ Structs to Flat Arrays)
    std::vector<double> h_pos_x(N), h_pos_y(N), h_pos_z(N), h_vel_x(N), h_vel_y(N), h_vel_z(N);

    for(int i=0; i<N; i++) {
        h_pos_x[i] = Particles[i].position.x;
        h_pos_y[i] = Particles[i].position.y;
        h_pos_z[i] = Particles[i].position.z;
        h_vel_x[i] = Particles[i].velocity.x;
        h_vel_y[i] = Particles[i].velocity.y;
        h_vel_z[i] = Particles[i].velocity.z;
    }

    size_t bytes = N * sizeof(double);
    
    // 2. Copy Data to GPU
    cudaMemcpy(mem.d_pos_x, h_pos_x.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(mem.d_pos_y, h_pos_y.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(mem.d_pos_z, h_pos_z.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(mem.d_vel_x, h_vel_x.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(mem.d_vel_y, h_vel_y.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(mem.d_vel_z, h_vel_z.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemset(mem.d_acc_x, 0, bytes); 
    cudaMemset(mem.d_acc_y, 0, bytes); 
    cudaMemset(mem.d_acc_z, 0, bytes);
    cudaMemset(mem.d_potential_energy, 0, bytes);
    cudaMemset(mem.d_kinetic_energy, 0, bytes);
}

// cudaFree all the memory
void free_gpu_memory(GpuMemory &mem) {
    cudaFree(mem.d_pos_x); cudaFree(mem.d_pos_y); cudaFree(mem.d_pos_z);
    cudaFree(mem.d_sorted_pos_x); cudaFree(mem.d_sorted_pos_y); cudaFree(mem.d_sorted_pos_z);
    cudaFree(mem.d_vel_x); cudaFree(mem.d_vel_y); cudaFree(mem.d_vel_z);
    cudaFree(mem.d_acc_x); cudaFree(mem.d_acc_y); cudaFree(mem.d_acc_z);
    cudaFree(mem.d_potential_energy);
    cudaFree(mem.d_kinetic_energy);

    cudaFree(mem.d_particle_indices);
    cudaFree(mem.d_cell_ids);
    cudaFree(mem.d_cell_start);
    cudaFree(mem.d_cell_end);
}

