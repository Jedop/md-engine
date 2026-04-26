#pragma once

#include <string>

struct SimConfig {
  int frames = 10000;
  double dt = 0.001;
  double rho = 1.0;
  bool fcc_or_not = true;
  int particles_per_side = 20;
  int unit_cells_per_side = 10;
  bool use_thermostat = true;
  double target_T = 1.5;
  bool do_annealing = false;
  bool do_time_reversal = false;
  std::string traj_file = "trajectory.xyz";
  std::string data_file = "data.txt";
};

void parse_args(int argc, char **argv, SimConfig &cfg);
void print_config(const SimConfig &cfg);
void print_help();
