#ifndef _FOCK_H_
#define _FOCK_H_

#include <complex.h>
#include "elementary.h"
#include <stdio.h>


/* Functions for working with single 'Fock' states, which 
 * are here defined as solutions to nointeracting Hamiltonian */
/* ========================================================== */

#define FOCK_CONFIG_STATE_CAPACITY 12 /* Max number of particles in fock state */

typedef struct {
	uint states[FOCK_CONFIG_STATE_CAPACITY]; /* Single particles states in sorted order 0, ..., n */
	uint size; /* Number of particles */
} Fock;

Fock fock_init(uint size);
Fock fock_init_direct(uint size, const uint* nrs);
Fock fock_init_quick(uint size, const int* m, const uint* rad_exc);

bool fock_equal(const Fock* s1, const Fock* s2);

double fock_energy(const Fock* fock);
int    fock_angular_momentum(const Fock* state);

uint   fock_operator_count(const Fock* fock, uint nr);
double fock_operator_create(Fock* fock, uint nr);
double fock_operator_annihilate(Fock* fock, uint nr);

complex double fock_compute(Fock* fock , double* x, double* y);
double fock_compute_density(Fock* fock , double* x, double* y);
double fock_compute_phase(Fock* fock , double* x, double* y);

double fock_matrix_element_kinetic(const Fock* left, const Fock* right);
double fock_matrix_element_delta_potential(const Fock* left, const Fock* right);


/* === Way of handling sets of fock states === */
/* =========================================== */

/* This structers allows indexing Fock states,
 * f.e. constructing linear combinations for same basis */

typedef struct {
	Fock* states;
	uint size;
        uint capacity;
} FockBasis;

FockBasis fock_basis_alloc(uint n_particles, double e_cutoff);
FockBasis fock_basis_alloc_fixed_angular_momentum(int m_tot, double e_cutoff, uint n_particles);
FockBasis fock_basis_alloc_copy(FockBasis src);

/* This functions DOES NOT check if fock is already in basis */
void fock_basis_add(FockBasis* basis, const Fock* fock);

/* Last state is moved to the place of removed state,
 * if idx is out of bounds, there is no error and
 * simply no state is removed */
void fock_basis_remove(FockBasis* basis, uint idx);

void fock_basis_free(FockBasis* basis);

/* Crude way of computing values of linear combinations of basis states,
 * probably better way is using FockState structere, you can create
 * such a structure from array of real numbers, see  (TODO: what is that function?) */

complex double fock_basis_compute(FockBasis set, double* x, double* y, const double* coeff);
double fock_basis_compute_density(FockBasis set, double* x, double* y, const double* coeff);
double fock_basis_compute_phase(FockBasis set, double* x, double* y, const double* coeff);

/* Returns set.size if element is not found */
uint fock_basis_find(FockBasis set, const Fock* state);

/* Those functions are mainly to see how your typical
 * perturbation procedure behaves for first few orders,
 * 2-nd and 3-rd order needs some limited basis */

double fock_perturbation_energy_one(double g, Fock *state);
double fock_perturbation_energy_two(double g, FockBasis set, Fock *state);
double fock_perturbation_energy_three(double g, FockBasis set, Fock *state);

typedef struct FockState {
        FockBasis basis;
        complex double* a;
} FockState;

void fock_state_normalize(FockState* fs);

/* Some coefficients 'a' in the FockState might be small,
 * you can keep those that |a[j]| >= cutoff, then size of basis
 * is also reduced (to safe space also capacity), the state is renormalized */
void fock_state_compress(FockState* state, double cutoff);

FockState fock_state_alloc(FockBasis set, double* coeff, double cutoff);
void fock_state_free(FockState* set);

complex double fock_state_compute(FockState cset, double* x, double* y);
complex double fock_state_aiaj(FockState set, uint i, uint j);

FockState fock_operator_psi(FockState fs, double x, double y);

typedef struct DEigenSol {
        complex double* U;
        double* val;
        uint size;
} DEigenSol;

typedef struct ZEigenSol {
        complex double* U;
        complex double* val;
        uint size;
} ZEigenSol;

void fock_matrix_fill(FockBasis set, double* T, double* V);
DEigenSol fock_eigen_vt_solve(double g, uint n, double* T, double* V);
DEigenSol fock_eigen_solve(double g, FockBasis set);




/* === Printing and saving to file === */
/* =================================== */

void state_fprint(FILE* f, uint state_nr);
static inline void fock_single_particle_print(uint state_nr) {
	state_fprint(stdout, state_nr);
}

void fock_fprint(FILE* f, const Fock* fock);
static inline void fock_print(const Fock* fock) {
	fock_fprint(stdout, fock);
}

void fock_basis_fprint(FILE* f, const FockBasis set);
static inline void fock_basis_print(const FockBasis set) {
	fock_basis_fprint(stdout, set);
}

void fock_basis_fprint_weighted(FILE* f, const FockBasis set, const double* coeff, uint n);
static inline void fock_basis_print_weighted(const FockBasis set, const double* coeff, uint n) {
	fock_basis_fprint_weighted(stdout, set, coeff, n);
}

void fock_coeffstateset_fprint(FILE* f, const FockState set, uint n);
static inline void fock_coeffstateset_print(const FockState set, uint n) {
	fock_coeffstateset_fprint(stdout, set, n);
}

/* === Some additional utility === */
/* =============================== */

double complex permanent_z(complex double* A, int n);

void util_bin_tofile(const void* src, uint nbytes, const char* filename);
void* util_bin_alloc_fromfile(uint nbytes, const char* filename);

#endif /* _FOCK_H_ */
