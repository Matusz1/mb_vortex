#ifndef _GPE_RADIAL_H_
#define _GPE_RADIAL_H_
#include "quick_function.h"
#include "fock.h"

#define BESSEL_LINEAR_COMB_MAX_SIZE 30

typedef unsigned int uint;

typedef struct {
        uint size;
        uint nrs[BESSEL_LINEAR_COMB_MAX_SIZE];
        double w[BESSEL_LINEAR_COMB_MAX_SIZE];
} BesselLinear;

typedef struct {
        bool isUp;
        uint nr;
} UVBesselState;

typedef struct {
        uint size;
        UVBesselState* types;
        double complex* evals;
        double complex* evecs;
        bool* plus_family;
} UVSol;

BesselLinear gpe_solve_m1(double gN, double* mi, double* Epp);
BesselLinear gpe_solve_m0(double gN, double* mi, double* Epp);
QuickFunction gpe_bes_lin_get_qf(const BesselLinear* bl);

UVSol bdg_solve(uint half_size, int m, double gN);
UVSol bdg_solve_m0(uint half_size, int m, double gN);
void uvsol_free(UVSol* uvsol);

#endif
