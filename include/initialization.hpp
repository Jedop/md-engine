#include "Particle.hpp"
#include "Vectors.hpp"
#include <cstdlib>
#include <ctime>
#include <vector>

std::vector<Particle> init_fcc_lattice(int unit_cells_per_side,
                                       double box_length);

std::vector<Particle> init_sc_lattice(int particles_per_side,
                                      double box_length);

void remove_center_of_mass_momentum(std::vector<Particle> &Particles);
