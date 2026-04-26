#pragma once

#include <fstream>
#include <iostream>
#include <omp.h>
#include <tuple>
#include <vector>

#include "Integrators.hpp"
#include "Particle.hpp"
#include "cell_list.hpp"
#include "config.hpp"
#include "forces.hpp"
#include "initialization.hpp"
#include "thermostat.hpp"

std::tuple<double, double, double> update(std::vector<Particle> &Particles,
                                          std::vector<int> &head,
                                          std::vector<int> &next, double dt,
                                          int nx, double cell_size, double box);

void run_simulation(SimConfig config);
