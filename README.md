# MD-Engine: A 3D C++ Molecular Dynamics Engine

## Overview

A 3D Molecular Dynamics Engine written in C++.

Currently, this project implements the following

- Create a 3D N body simulation framework in C++ with Periodic Boundary Conditions(PBC).
- Implement a standard Lennard-Jones (LJ) potential to model a simple nobel gas.
- Implement multiple numerical integrators, specifically: Forward Euler, Velocity Verlet, and DKD Leapfrog methods.

## Key Features

- Cell lists: Divides the space into cubes to optimize the Force Calculation Algorithm from $O(N^2)$ to $O(N)$, enabling simulation of 8000 particles for 10,000 timesteps in 170s.
- Multiple Integrators: This enables comparison and numerical analyses of various methods of integration.

## Physics & Implementation Details

To ensure physical accuracy and numerical stability, the engine implements standard molecular dynamics conventions:

- **Lennard-Jones Reduced Units:** All calculations (mass, energy, distance, time) are performed in dimensionless LJ reduced units to prevent floating-point underflow/overflow.
- **Initialization:** Particles are initialized in a stable 3D lattice configuration. Initial velocities are assigned randomly centered at 0 to avoid large drift velocities.
- **Minimum Image Convention (MIC):** Implemented for calculating the shortest distance between particles across Periodic Boundary Conditions without expensive division operations.
- **Force Truncation:** The LJ potential is cut off at $r_c = 2.5\sigma$ to optimize calculations.

## Numerical Analysis

### Euler vs Leapfrog vs Velocity Verlet

| Integrator       | Euler                            | Leapfrog                               | Velocity Verlet                                     |
| ---------------- | -------------------------------- | -------------------------------------- | --------------------------------------------------- |
| Energy over time | ![Euler](assets/EulerEnergy.png) | ![Leapfrog](assets/LeapfrogEnergy.png) | ![Velocity Verlet](assets/VelocityVerletEnergy.png) |

- Euler : Energy increases exponentially over time, since it is **not** a symplectic integrator. This rules out Euler integration for this project.
- Leapfrog: Energy is conserved over time, since it is a symplectic integrator.
- Velocity Verlet: Energy is also conserved in this case, as it is a symplectic integrator as well.

#### Why Velocity Verlet over Leapfrog?

Clearly, the difference in energy drift between Leapfrog and Velocity Verlet Integration is close to negligible. One might wonder then, why choose one over the other?

In our case, it is very clear that Velocity Verlet is the superior choice. This is because the Leapfrog method is such that at each timestep n, the calculated positions(x) are at timestep n, but the calculated velocities(v) are at timestep $(n + 1/2)$. This means, every time we want to calculate the kinetic energy(or any quantity involving the velocities of the particles), we need to take it back a half-step outside of the main loop, and **then** calculate the required quantities. This is clearly suboptimal. This is a non-issue for Velocity Verlet, as it always works with a full step of velocity and position at the end of each timestep. This makes it the superior choice for Molecular Dynamics, where quantities involving position and velocity must be calculated often.

## Visualization

The engine outputs raw trajectory data in the standard `.xyz` file format. This allows for seamless integration with scientific visualization software like **OVITO**.

## ![Output](assets/Visualization.gif)
