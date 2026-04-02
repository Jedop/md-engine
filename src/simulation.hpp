#pragma once

#include <tuple>
#include <vector>

#include "Integrators.hpp"
#include "Particle.hpp"
#include "cell_list.hpp"
#include "forces.hpp"

std::tuple<double, double, double> update(std::vector<Particle> &Particles,
                                          std::vector<int> &head,
                                          std::vector<int> &next, double dt);
