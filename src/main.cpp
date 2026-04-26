#include "simulation.hpp"

int main() {

  SimConfig cfg = {20000, 0.001, 1, true, 20, 10, true, 1.5, true, false};

  run_simulation(cfg);

  return 0;
}
