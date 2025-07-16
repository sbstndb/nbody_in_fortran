# N-Body Simulation in Fortran

This repository contains a simple N-body gravitational simulation code written in Fortran. It simulates the motion of a system of particles under their mutual gravitational attraction.

## Code Structure

The code is organized into two main files:

-   `module.f90`: This file defines the `mod_compute` module, which contains the core logic for the simulation.
    -   `particle_s`: A derived data type to store the position (`x`, `y`, `z`) and velocity (`vx`, `vy`, `vz`) of each particle.
    -   `init`: A subroutine to initialize the particles' positions and velocities to zero.
    -   `move`: The main integration subroutine. It calculates the gravitational forces between all pairs of particles and updates their positions and velocities.
    -   `loop`: A subroutine that calls `move` for a specified number of time steps.
    -   `compute`: The main driver subroutine that allocates memory for particles, initializes them, runs the simulation loop, and prints the final position of the first particle.
-   `main.f90`: The main program that calls the `compute` subroutine from the `mod_compute` module.

## Equations of Motion

The simulation uses a direct, brute-force N-body algorithm, where the force on each particle is computed by summing the gravitational influences of all other particles.

The gravitational force \(\vec{F}_{ij}\) exerted by a particle \(j\) on a particle \(i\) is given by Newton's law of universal gravitation:

```math
\vec{F}_{ij} = G \frac{m_i m_j}{|\vec{r}_{ij}|^2} \frac{\vec{r}_{ij}}{|\vec{r}_{ij}|}
```

where:
- \(G\) is the gravitational constant.
- \(m_i\) and \(m_j\) are the masses of particles \(i\) and \(j\).
- \(\vec{r}_{ij} = \vec{r}_j - \vec{r}_i\) is the displacement vector from particle \(i\) to particle \(j\).

In this specific implementation, the masses (\(m_i\), \(m_j\)) and the gravitational constant (\(G\)) are assumed to be 1 for simplicity.

The total force on particle \(i\) is the vector sum of the forces from all other particles:

```math
\vec{F}_i = \sum_{j=1, j \neq i}^{N} \vec{F}_{ij}
```

The acceleration of particle \(i\) is then \(\vec{a}_i = \vec{F}_i / m_i\).

To avoid numerical instability when two particles get very close, a **softening factor** (\(\epsilon\)) is introduced to the distance calculation. The denominator becomes:

```math
(|\vec{r}_{ij}|^2 + \epsilon^2)^{3/2}
```

This prevents the force from becoming infinite. In the code (`subroutine move`), this is implemented as `d_2 = dx_t*dx_t + dy_t*dy_t + dz_t*dz_t + softening`.

The time integration is performed using a **Leapfrog (kick-drift-kick)** scheme, which is a second-order method common for N-body simulations. In `subroutine move`:
1.  **Kick**: Velocities are updated based on the forces at the current positions.
    ```fortran
    part(i)%vx = part(i)%vx + dt*fx_t
    ```
2.  **Drift**: Positions are updated using the new velocities.
    ```fortran
    part(i)%x = part(i)%x + dt*part(i)%vx
    ```

## Compilation

The code is compiled using the Intel Fortran Compiler (`ifort`). A simple build script is provided.

To compile the code, run:
```bash
sh make.sh
```
This will generate an executable named `main`.

To run the simulation:
```bash
./main
```

## TODO

Here are three suggestions for improving this code:

- [ ] **Parallelization**: The main computational bottleneck is the nested loop in the `move` subroutine, which has a complexity of \(O(N^2)\). This loop is highly parallelizable.
    - [ ] **OpenMP**: Introduce OpenMP directives to parallelize the outer loop over particles (`do i = 1, n`). This would be a straightforward shared-memory parallelization strategy. Care must be taken with the reduction of forces.
    - [ ] **MPI**: For larger-scale simulations distributed across multiple nodes, MPI could be used to decompose the particle data among processes.

- [ ] **Algorithmic Optimization**: The \(O(N^2)\) complexity is inefficient for a large number of particles.
    - [ ] **Tree-based methods**: Implement a Barnes-Hut or Fast Multipole Method (FMM) to reduce the complexity to \(O(N \log N)\) or \(O(N)\) respectively. These methods approximate the force from distant clusters of particles instead of computing direct interactions.
    - [ ] **Neighbor Lists**: For systems with short-range interactions, Verlet lists or cell lists can be used to only compute forces for nearby particles.

- [ ] **IO and Visualization**: The current output is minimal.
    - [ ] **Data Output**: Implement proper data output using standard file formats like HDF5 or NetCDF. This would allow for storing the state of the system at different time steps for analysis.
    - [ ] **Visualization**: Connect the output to a visualization tool like ParaView or Visit. Writing data in a format they can read (like VTK) would enable visual analysis of the particle trajectories.
    - [ ] **Input Parameters**: Allow setting simulation parameters (number of particles, time steps, softening factor) via a command-line argument or an input file instead of hardcoding them.
