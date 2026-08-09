#include "config.hpp"
#include "simulation.hpp"

int main(int argc, char **argv) {
  SimConfig cfg;

  parse_args(argc, argv, cfg);

  run_simulation(cfg);
}
