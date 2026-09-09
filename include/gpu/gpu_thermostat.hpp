#pragma once
#include <vector>
#include "Particle.hpp"
#include "gpu_memory.hpp"
#include <iostream>

void apply_berendsen_thermostat_gpu(std::vector<Particle> &Particles,
                                double T_current, double T_target, double dt, GpuMemory &mem);

void reverse_velocities_gpu(int N, GpuMemory& mem);