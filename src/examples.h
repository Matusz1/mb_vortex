#ifndef _EXAMPLES_H_
#define _EXAMPLES_H_
#include "elementary.h"

void example_random_sampling(double gamma, uint npart);
void example_gpe_energy(uint npart);
void example_bdg_write(int m, double gN, bool vortex_state);
void example_bdg(uint npart, double gamma_max, uint ncuts);
void example_rho_diag(uint npart, double gamma_min, double gamma_max, uint ncuts, uint ssize);
void example_rho_evec(uint npart, double gamma, uint ssize);
void example_rho_gpe_compariton(uint npart, double gamma, uint ssize);
void example_bdg_uv(uint npart, double gamma);

#endif // !_EXAMPLES_H_
