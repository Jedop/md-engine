#include "config.hpp"
#include <iostream>

void print_help() {
  std::cout << "Usage: ./MD_Engine [options]\n\n"
            << "--dt <value>\n"
            << "--frames <value>\n"
            << "--rho <value>\n"
            << "--fcc / --sc\n"
            << "--particles-per-side <int>\n"
            << "--unit-cells <int>\n"
            << "--turn-off-thermostat\n"
            << "--target-T <value>\n"
            << "--anneal\n"
            << "--time-reversal\n"
            << "--traj <file>\n"
            << "--data <file>\n";
}

void print_config(const SimConfig &cfg) {
  std::cout << "Running with config:\n";
  std::cout << "frames=" << cfg.frames << ", dt=" << cfg.dt
            << ", rho=" << cfg.rho
            << ", lattice=" << (cfg.fcc_or_not ? "FCC" : "SC")
            << ", particles_per_side=" << cfg.particles_per_side
            << ", unit_cells=" << cfg.unit_cells_per_side
            << ", turn-off-thermostat=" << cfg.use_thermostat
            << ", target_T=" << cfg.target_T
            << ", annealing=" << cfg.do_annealing
            << ", time_reversal=" << cfg.do_time_reversal
            << ", traj=" << cfg.traj_file << ", data=" << cfg.data_file
            << "\n\n";
}

void parse_args(int argc, char **argv, SimConfig &cfg) {

  if (argc == 1) {
    print_config(cfg); // ← print defaults
    return;
  }

  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];

    if (arg == "-h" || arg == "--help") {
      print_help();
      exit(0);
    }

    else if (arg == "--dt" && i + 1 < argc)
      cfg.dt = std::stod(argv[++i]);

    else if (arg == "--frames" && i + 1 < argc)
      cfg.frames = std::stoi(argv[++i]);

    else if (arg == "--rho" && i + 1 < argc)
      cfg.rho = std::stod(argv[++i]);

    else if (arg == "--fcc")
      cfg.fcc_or_not = true;

    else if (arg == "--sc")
      cfg.fcc_or_not = false;

    else if (arg == "--particles-per-side" && i + 1 < argc)
      cfg.particles_per_side = std::stoi(argv[++i]);

    else if (arg == "--unit-cells" && i + 1 < argc)
      cfg.unit_cells_per_side = std::stoi(argv[++i]);

    else if (arg == "--turn-off-thermostat")
      cfg.use_thermostat = false;

    else if (arg == "--target-T" && i + 1 < argc) {
      cfg.target_T = std::stod(argv[++i]);
      cfg.use_thermostat = true;
    } else if (arg == "--anneal") {
      cfg.do_annealing = true;
      cfg.use_thermostat = true;
    } else if (arg == "--time-reversal")
      cfg.do_time_reversal = true;

    else if (arg == "--traj" && i + 1 < argc)
      cfg.traj_file = argv[++i];

    else if (arg == "--data" && i + 1 < argc)
      cfg.data_file = argv[++i];
  }

  print_config(cfg);
}
