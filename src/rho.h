#ifndef _RHO_H_
#define _RHO_H_

#include "fock.h"
#include "state.h"

typedef struct Rho {
        Basis basis;
        double complex* M;
} Rho;

Rho rho_alloc(FockState fs);
void rho_free(Rho* rho);

complex double rho_compute(Rho rho, double r, double phi, double rp, double phip);
complex double rho_compute_xy(Rho rho, double x, double y, double xp, double yp);
double rho_compute_diag(Rho rho, double r, double phi);
double rho_compute_diag_xy(Rho rho, double x, double y);

void deigen_sol_free(DEigenSol* es);
void zeigen_sol_free(DEigenSol* es);

DEigenSol rho_diag(Rho* rho);

/* TODO: How to compute diagonalized rho elements */

double* fock_state_alloc_random(FockState set, uint n);

typedef enum {
        FOCK_RHO_RAND_KEEP_ARR,
        FOCK_RHO_RAND_RELOAD,
} RhoRandOpt;

void rho_rand(Rho rho, RhoRandOpt opt, double* x, double *y);

#endif // !_RHO_H_
