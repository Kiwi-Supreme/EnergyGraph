# EnergyGraph

This C++ program simulates a simple molecular dynamics system of particles in a cubic lattice using the Lennard-Jones potential. The simulation includes initializing particle positions and velocities, calculating forces, and integrating particle movements over time. The program outputs the potential, kinetic, and total energies of the system to a file named `energy_data.txt`. Here is a breakdown of the key components:

1. **Initialization**: Particles are initialized in a cubic lattice structure with random velocities, which are adjusted to achieve the desired temperature.
2. **Force Calculation**: The forces between particles are computed using the Lennard-Jones potential, considering periodic boundary conditions.
3. **Integration**: The particle positions and velocities are updated using the velocity-Verlet algorithm.
4. **Energy Calculation**: The kinetic, potential, and total energies are computed at each time step.
5. **File Output**: The energy data is written to `energy_data.txt` at each time step.

### Key Functions
- `init()`: Initializes particle positions and velocities.
- `Force()`: Calculates forces between particles and updates potential energy.
- `integrate(double en)`: Updates particle positions and velocities, and computes kinetic and total energies.
- `ran2(long *idum)`: Generates random numbers for initializing velocities.

### Main Loop
- Runs the simulation for a specified maximum time (`tmax`) with a defined time step (`dt`).
- Computes forces and integrates particle positions at each step.
- Writes energy data to `energy_data.txt` for later analysis and visualization with tools like `gnuplot`.
