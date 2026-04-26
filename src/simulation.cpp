#include "simulation.hpp"
#include "utils.hpp"

std::tuple<double, double, double> update(std::vector<Particle> &Particles,
                                          std::vector<int> &head,
                                          std::vector<int> &next, double dt,
                                          int nx, double cell_size,
                                          double box) {
  // 1. Pre-Force Update (Move everyone to new positions)
  VelocityVerlet::step1(Particles, dt);

  // 2. Rebuild the Grid
  build_cell_lists(Particles, head, next, nx, cell_size);

  // 3. Synchronous Force Calculation (Calculate forces at the NEW positions)
  auto [new_acc, potential_energy] =
      compute_all_forces(Particles, head, next, nx, cell_size, box);

  // 4. Post-Force Update (Finish the time step using the new forces)
  VelocityVerlet::step2(Particles, new_acc, dt);

  // Calculates Kinetic Energy
  double kinetic_energy = 0;

  for (const auto &p : Particles) {
    Vec3 v = p.velocity;
    // v = p.velocity - p.acceleration * dt * 0.5; Leapfrog Calculation
    double v2 = v.x * v.x + v.y * v.y + v.z * v.z;
    kinetic_energy += 0.5 * v2;
  }

  double Energy = potential_energy + kinetic_energy;
  return {potential_energy, kinetic_energy, Energy};
}

void run_simulation(SimConfig config) {
  auto start = std::chrono::steady_clock::now();
  if (config.use_thermostat && config.do_time_reversal) {
    std::cout
        << "Don't try to use thermostat and do time reversal together, do "
           "either one. For now, it'll act as if the thermostat is off.";
  }

  // Constants
  const int N = config.fcc_or_not ? 4 * std::pow(config.unit_cells_per_side, 3)
                                  : std::pow(config.particles_per_side, 3);
  const double box = pow(N / config.rho, 1.0 / 3.0);
  const int nx = int(box / rc);
  const int num_cells = nx * nx * nx;
  const double cell_size = box / nx;

  // Initialization
  std::vector<Particle> Particles;

  if (config.fcc_or_not) {
    Particles = init_fcc_lattice(config.unit_cells_per_side, box);
  } else {
    Particles = init_sc_lattice(config.particles_per_side, box);
  }

  // This function also adds random velocities but I didn't know how to fit all
  // that information into one function name
  remove_center_of_mass_momentum(Particles);

  // head and next vectors act as linked list
  std::vector<int> head(num_cells, -1);
  std::vector<int> next(Particles.size());

  // Files into which everything is written
  std::ofstream traj_file("trajectory.xyz");
  std::ofstream data_file("thermo.dat");

  // Initial build of cell list and other initial stuff
  build_cell_lists(Particles, head, next, nx, cell_size);
  auto [initial_acc, _] =
      compute_all_forces(Particles, head, next, nx, cell_size, box);

  for (size_t i = 0; i < Particles.size(); i++) {

    Particles[i].acceleration = initial_acc[i];
    // Leapfrog Initialization
    // Particles[i].velocity += Particles[i].acceleration * dt * 0.5;
  }

  // Recording initial positions in case of time reversal
  std::vector<Vec3> initial_positions(Particles.size());
  for (size_t i = 0; i < Particles.size(); i++) {
    initial_positions[i] = Particles[i].position;
  }

  // Main Loop
  for (int step = 0; step < config.frames; step++) {

    // --- FEATURE: TIME REVERSAL ---
    if (config.do_time_reversal && step == (config.frames / 2)) {
      for (auto &p : Particles)
        p.velocity = p.velocity * -1.0;
      std::cout << ">>> TIME REVERSAL TRIGGERED <<<\n";
    }

    // --- FEATURE: ANNEALING ---
    if (config.do_annealing) {
      if (step < config.frames * 0.75) { // Cool down over 75% of the run
        config.target_T = 1.5 - (1.4 * ((double)step / (config.frames * 0.75)));
      } else {
        config.target_T = 0.1;
      }
    }

    // --- CORE PHYSICS ---
    auto [pot_E, kin_E, Tot_E] =
        update(Particles, head, next, config.dt, nx, cell_size, box);

    // --- FEATURE: THERMOSTAT ---
    if (config.use_thermostat && !config.do_time_reversal) {
      // Safety: Never use thermostat during time reversal!
      double T_current = (2 * kin_E) / (3 * Particles.size());
      apply_berendsen_thermostat(Particles, T_current, config.target_T,
                                 config.dt);
    }

    if (step % 100 == 0 && step > 0) {

      auto now = std::chrono::steady_clock::now();
      double elapsed = std::chrono::duration<double>(now - start).count();
      double progress = (double)step / config.frames;
      double eta = elapsed * (1.0 / progress - 1.0);
      print_progress(step, config.frames, eta);

      traj_file << Particles.size() << "\n";

      traj_file << "Lattice=\"" << box << " 0.0 0.0 0.0 " << box
                << " 0.0 0.0 0.0 " << box
                << "\" Properties=species:S:1:pos:R:3\n";

      for (const auto &p : Particles) {
        traj_file << "Ar " << p.position.x << " " << p.position.y << " "
                  << p.position.z << "\n";
      }

      data_file << pot_E << " " << kin_E << " " << Tot_E << " "
                << (2 * kin_E / (3 * N)) << "\n";
    }
  }
  traj_file.close();
  data_file.close();

  print_progress(config.frames, config.frames, 0);
  std::cout << "\n";

  // Output Time Reversal Error
  if (config.do_time_reversal) {
    double max_error = 0.0;

    for (size_t i = 0; i < Particles.size(); i++) {
      Vec3 diff = Particles[i].position - initial_positions[i];

      // PBC is unnecessary because ideally they'll be at initial positions
      double error_distance =
          std::sqrt(diff.x * diff.x + diff.y * diff.y + diff.z * diff.z);

      if (error_distance > max_error) {
        max_error = error_distance;
      }
    }

    std::cout << "====================================\n";
    std::cout << "Time Reversibility Maximum Absolute Error: " << max_error
              << "\n";
    std::cout << "====================================\n";
  }
}
