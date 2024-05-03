#ifndef _STATE_H_
#define _STATE_H_

#include "elementary.h"
#include "quick_function.h"
#include <complex.h>

/* Single particle bessel states, they are handled using
 * uint's, and for more detail about a state you can call
 * appropriate function (state_info(indx)) */

double state_energy(uint n);

double state_value_r(uint n, double r);
complex double state_value(uint n, double r, double phi);
complex double state_value_xy(uint n, double x, double y);

/* Ensures OpenMP does not break when evaulating
 * values of states up to, including max */
void state_update_quick_functions(uint max);

/* Info about state, only retured from functions */

typedef struct {
	int m;
	uint zero_num;
        double zero_val;
        double normalization_constant;
} StateInfo;

StateInfo state_info(uint n);
StateInfo state_info_direct(int m, uint zero_num);
uint state_find_index(StateInfo state);

const QuickFunction* state_get_qf(uint nr);



/* When we want to work with linear combinations of single particle states */

typedef struct {
        uint* states;
        uint size;
        uint capacity;
} Basis;

void state_basis_free(Basis* basis);
bool state_basis_contains(Basis* basis, uint state);
Basis state_basis_alloc_copy(Basis* basis);

/* Does not check whether the state is already in */
void state_basis_add(Basis* basis, uint state);

#endif // _STATE_H_
