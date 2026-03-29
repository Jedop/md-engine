#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

constexpr double box = 70;
constexpr double box_r = 1 / box;

struct Vec3 {
  double x, y, z;

  Vec3 operator+(const Vec3 &o) const { return {x + o.x, y + o.y, z + o.z}; }

  Vec3 operator-(const Vec3 &o) const { return {x - o.x, y - o.y, z - o.z}; }

  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }

  Vec3 &operator+=(const Vec3 &o) {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }

  void apply_pbc(double box, double box_r) {
    auto wrap = [&](double &d) {
      int k = (int)(d * box_r + (d >= 0 ? 0.5 : -0.5));
      d -= k * box;
    };
    wrap(x);
    wrap(y);
    wrap(z);
  }
};

struct Particle {
  Vec3 position;
  Vec3 velocity;
  Vec3 acceleration;

  Particle(double x, double y, double z) {
    position = {x, y, z};
    velocity = {0, 0, 0};
    acceleration = {0, 0, 0};
  }
};

std::vector<Vec3> compute_all_forces(const std::vector<Particle> &Particles) {
  std::vector<Vec3> all_acc(Particles.size(), {0, 0, 0});

  for (size_t i = 0; i < Particles.size(); i++) {
    for (size_t j = 0; j < Particles.size(); j++) {
      if (i == j)
        continue; // Skip self

      Vec3 r = Particles[i].position - Particles[j].position;
      r.apply_pbc(box, box_r);
      double r2 = r.x * r.x + r.y * r.y + r.z * r.z; // r^2
      double rc = 2.5 * 0.01;

      if (r2 < 1e-12 || r2 > rc * rc)
        continue; // Safety

      double inv_r2 = 1.0 / r2;
      double sigma2 = 0.01; // sigma^2
      double inv_r6 = pow(sigma2 * inv_r2, 3);

      all_acc[i] += r * 24 * 0.1 * inv_r2 * inv_r6 * (2 * inv_r6 - 1);
    }
  }
  return all_acc;
}
namespace VelocityVerlet {

void step1(std::vector<Particle> &Particles, double dt) {
  for (auto &p : Particles) {
    p.position += p.velocity * dt + p.acceleration * (dt * dt / 2.0);
    p.velocity += p.acceleration * (dt / 2.0);
  }
}

void step2(std::vector<Particle> &Particles,
           const std::vector<Vec3> &new_accelerations, double dt) {
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].acceleration = new_accelerations[i];
    Particles[i].velocity += Particles[i].acceleration * (dt / 2.0);
  }
}

} // namespace VelocityVerlet

void update(std::vector<Particle> &Particles, double dt) {
  // 1. Pre-Force Update (Move everyone to new positions)
  VelocityVerlet::step1(Particles, dt);

  // 2. Synchronous Force Calculation (Calculate forces at the NEW positions)
  std::vector<Vec3> new_acc = compute_all_forces(Particles);

  // 3. Post-Force Update (Finish the time step using the new forces)
  VelocityVerlet::step2(Particles, new_acc, dt);
}

int main() {
  std::vector<Particle> Particles;
  Particles.reserve(512);

  int particles_per_side = 8;
  double spacing = box / particles_per_side;

  for (int i = 0; i < particles_per_side; i++) {
    for (int j = 0; j < particles_per_side; j++) {
      for (int k = 0; k < particles_per_side; k++) {

        double x = (i + 0.5) * spacing;
        double y = (j + 0.5) * spacing;
        double z = (k + 0.5) * spacing;

        Particles.emplace_back(x, y, z);
      }
    }
  }
  Particles[3].velocity = {1.0, 1.0, 1.0};
  std::vector<Vec3> initial_acc = compute_all_forces(Particles);
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].acceleration = initial_acc[i];
  }

  int no_of_frames = 20000;
  std::ofstream traj_file("trajectory.xyz");

  for (int i = 0; i < no_of_frames; i++) {
    std::cout << i << "\n";
    update(Particles, 0.001);
    if (i % 1 == 0) {
      traj_file << Particles.size() << "\n";

      // THE MAGIC LINE: This tells OVITO you have a cubic box of size L x L x L
      traj_file << "Lattice=\"" << box << " 0.0 0.0 0.0 " << box
                << " 0.0 0.0 0.0 " << box
                << "\" Properties=species:S:1:pos:R:3\n";

      for (const auto &p : Particles) {
        // Output the actual, UNWRAPPED drifting coordinates
        traj_file << "Ar " << p.position.x << " " << p.position.y << " "
                  << p.position.z << "\n";
      }
    }
  }
  traj_file.close();
  return 0;
}
