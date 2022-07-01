# Simulation of a 1D Quantum Harmonic Oscillator using the Metropolis Hasting Algorithm
## Authors: Nicolás A. N. Salas, Fabian Felipe Quevedo, Julian Mateo de Mendoza, Jefferson Garzón
After cloning the repository, it is enough to use make in the main folder to compile the proyect, it will produce an executable named foo. The executable receives as input:
- time steps
- Metropolis Hasting steps
- mass
to use the parallelization it is necessary to run it using mpirun. An example with 4 cores is:
mpirun -np 4 ./foo 1000 5000 1.

The tests can be run using make tests, and you can also clean all the built files using make clean.

