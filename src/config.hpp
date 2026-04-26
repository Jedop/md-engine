#pragma once

struct SimConfig {
  int frames;
  double dt;
  double rho;
  bool fcc_or_not;
  int particles_per_side;
  int unit_cells_per_side;
  bool use_thermostat;
  double target_T;
  bool do_annealing;
  bool do_time_reversal;
};
