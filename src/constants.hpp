#pragma once

constexpr int no_of_frames = 10000;
constexpr double dt = 0.001;
constexpr int particles_per_side = 20;
constexpr int N = particles_per_side * particles_per_side * particles_per_side;
constexpr double box = particles_per_side * 1.122;
constexpr double box_r = 1 / box;
constexpr double rc = 2.5;

constexpr int nx = int(box / rc);
constexpr int num_cells = nx * nx * nx;
constexpr double cell_size = box / nx;
