# Solver for many-body boson systems (WIP)

Software for exact diagonalization of contact interaction bosons in 2D cylinder,
density matrix diagonalization, positions drawing, self-consistent solving special case of 
GPE for centered vortex in cylinder. Currenty working on simple BdG solver fot vortex state from GPE.

* Many body calculations are being done with Fock state basis consisting of single particle states that are solutions of single particle Schrödinger equation (Bessel functions).
* Many body function is used to determine density matrix (this we can compare to mean-field solutions).
* Density matrix can be used to draw positions of particles one-by-one.
* Solving GPE equation with centered vortex can be done very precisly using 1st bessel function and its zeros.
* Solution to GPE from previous points can be used in BdG equation (WIP).

For bit more details contact publication: ... (WIP).
