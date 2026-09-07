#pragma once

#include <cmath>
#include <omp.h>
#include <vector>

#include "Particle.hpp"
#include "Vectors.hpp"
#include "cell_list.hpp"
#include "constants.hpp"

void apply_velocity_rescaling(std::vector<Particle> &Particles,
                              double T_current, double T_target, double dt);

void apply_berendsen_thermostat(std::vector<Particle> &Particles,
                                double T_current, double T_target, double dt);
