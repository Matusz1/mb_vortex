#ifndef _FOCK_H_
#define _FOCK_H_

#include <stdbool.h>
#include <complex.h>

/* For option for printing out states */
#include <stdio.h>

typedef unsigned int uint;

#define FOCK_CONFIG_STATE_CAPACITY 6 /* Max number of particles in fock state */
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

typedef enum {
	STATE_PLANEWAVE,
	STATE_BESSEL,
	_STATE_NUM
} StateType;

typedef struct {
	uint states[FOCK_CONFIG_STATE_CAPACITY];   /* Single particles states in sorted order 0, ..., n */
	uint size;      /* Number of particles */
	StateType type; /* Type of state fe. plane waves */
} FockState;

typedef struct {
	int kx;
	int ky;
} Momentum;

typedef struct Vec2z {
        complex double x;
        complex double y;
} Vec2z;

typedef struct Vec2d {
        double x;
        double y;
} Vec2d;

typedef struct Lattice2D {
        double* x;
        double* y;
        uint size;
} Lattice2D;

/* Do not initialize this structure by yourself,
 * use funtion 'fock_bessl_wave_state(int m, uint zero_num)' */
typedef struct {
	int m;
	uint zero_num;
        double zero_val; /* Automaticaly created when the state is created */
        double normalization_constant; /* Automaticaly created when the state is created */
} BesselWaveState;

typedef struct {
        uint* is;
        uint* js;
        complex double* coeff;
        uint size;
        StateType type;
} G1State;

typedef Momentum PlaneWaveState;

#define MOMENTUM_INIT(kx, ky) { kx, ky }
#define PLANE_WAVE_STATE_INIT(kx, ky) { kx, ky }

double complex permanent_z(complex double* A, int n);

Lattice2D fock_lattice_alloc(Vec2d ll, Vec2d ur, uint n_rows, uint n_cols);
Lattice2D fock_lattice_alloc_helix(uint n_points, double r);
void fock_lattice_free(Lattice2D* lattice);

/* Some general functions for fock states */
/* ====================================== */

FockState fock_init(uint size, StateType type);
FockState fock_init_direct(uint size, StateType type, const uint* nrs);
FockState fock_init_quick_bessel(uint size, const int* m, const uint* rad_exc);

bool fock_equal(const FockState* s1, const FockState* s2);

double fock_energy(const FockState* fock);

uint fock_operator_count(const FockState* fock, uint nr);
double fock_operator_create(FockState* fock, uint nr);
double fock_operator_annihilate(FockState* fock, uint nr);

/* Functions for particular type of single particle fucntions */
/* ========================================================== */

double fock_single_energy(uint state_nr, StateType type);
complex double fock_single_compute(uint state_nr, StateType type, double x, double y);

/* Plane wave states */

PlaneWaveState fock_planewave_state_index(uint state_nr);
double fock_planewave_energy(uint state_nr);
uint fock_planewave_find_index(PlaneWaveState state);
Momentum fock_planewave_momentum(const FockState* fock);

/* States in cylinder - Bessel functions */

BesselWaveState fock_bessel_state_index(uint state_nr);
BesselWaveState fock_bessel_state_direct(int m, uint zero_num);
uint fock_bessel_find_index(BesselWaveState state);
int fock_bessel_angular_momentum(const FockState* state);
double bessel_value_radial(uint state_nr, double r);
double complex bessel_value(uint state_nr, double r, double phi);

double fock_matrix_element_delta_potential(const FockState* left, const FockState* right);
double fock_matrix_element_kinetic(const FockState* left, const FockState* right);


/* === Way of handling sets of fock states === */
/* =========================================== */

typedef struct {
	FockState* states;
	uint size;
} FockStateSet;

typedef struct FockCoeffStateSet {
        FockStateSet set;
        complex double* coeff;
} FockCoeffStateSet;


FockStateSet fock_stateset_alloc(StateType type, uint n_states, uint n_particles, double e_cutoff);
FockStateSet fock_stateset_alloc_same_momentum(int kx, int ky, int e_cutoff, uint max_state, uint n_particles);
FockStateSet fock_stateset_alloc_same_angular_momentum(int m_tot, double e_cutoff, uint max_index, uint n_particles);
FockStateSet fock_stateset_alloc_copy(FockStateSet src);
void fock_stateset_append(FockStateSet* set, const FockState* fock);
void fock_stateset_free(FockStateSet* set);

double* fock_stateset_perturbation_function(FockStateSet left, FockState* state);
double fock_stateset_perturbation_energy_one(double g, FockState *state);
double fock_stateset_perturbation_energy_two(double g, FockStateSet set, FockState *state);
double fock_stateset_perturbation_energy_three(double g, FockStateSet set, FockState *state);

FockCoeffStateSet fock_operator_psi(FockCoeffStateSet cset, double x, double y);
FockCoeffStateSet fock_coeffstateset_alloc_dcoeff(FockStateSet set, double* coeff, double cutoff);
void fock_coeffstateset_free(FockCoeffStateSet* set);
complex double fock_coeffstateset_compute(FockCoeffStateSet cset, double* x, double* y);
complex double fock_coeffstateset_aiaj(FockCoeffStateSet set, uint i, uint j);

/* Returns set.size if element is not found */
uint fock_stateset_find(FockStateSet set, const FockState* state);

typedef struct EigenSolutionSet {
	double* values;
	double* vectors;
	uint size;
} EigenSolutionSet;

typedef struct CmplxEigenSolutionSet {
	double* values;
	complex double* vectors;
	uint size;
} CmplxEigenSolutionSet;

complex double fock_compute(FockState* fock , double* x, double* y);
Vec2z fock_compute_grad(FockState* fock, double* x, double* y, uint n);
double fock_compute_density(FockState* fock , double* x, double* y);
double fock_compute_phase(FockState* fock , double* x, double* y);

complex double fock_stateset_compute_weighted(FockStateSet set, double* x, double* y, const double* coeff);
Vec2z fock_stateset_compute_weighted_grad(FockStateSet set, double* x, double* y, const double* coeff, uint n);
double fock_stateset_compute_density_weighted(FockStateSet set, double* x, double* y, const double* coeff);
double fock_stateset_compute_phase_weighted(FockStateSet set, double* x, double* y, const double* coeff);
void fock_stateset_bohm_next(FockStateSet set, double* x, double* y, const double* coeff, double dt);

double* fock_stateset_metropolis_alloc(FockStateSet set, const double* coeff, uint n_vals, uint n_bournout, unsigned long seed);

double fock_stateset_aiaj(FockStateSet set, const double* coeff, uint i, uint j);
G1State fock_g1_alloc(FockStateSet set, const double* coeff);
G1State fock_g1_alloc_cset(FockCoeffStateSet set);
double* fock_coeffstateset_alloc_random(FockCoeffStateSet set, uint n);
void fock_g1_free(G1State* g1);
double complex fock_g1_compute(G1State* g1, double x, double y, double xp, double yp);

typedef enum {
        FOCK_G1_RAND_KEEP_ARR,
        FOCK_G1_RAND_RELOAD,
} FockG1RandOpt;

Vec2d fock_g1_rand(G1State* g1, FockG1RandOpt opt);

complex double* fock_g1_matrix_alloc(Lattice2D lattice, G1State* g1);
double* fock_g1_matrix_calc_top_eigenvector(double complex* M, uint n, uint n_iter, double* eval);
CmplxEigenSolutionSet fock_g1_matrix_solve(double complex* M, uint dim);
void fock_eigen_save(const char* fname, EigenSolutionSet set);
EigenSolutionSet fock_eigen_alloc_file(const char* fname);

void fock_matrix_fill(FockStateSet set, double* T, double* V);
EigenSolutionSet fock_eigen_vt_solve(double g, uint n, double* T, double* V);
EigenSolutionSet fock_eigen_omega_vt_solve(double g, double omega, FockStateSet fset, double* T, double* V);
EigenSolutionSet fock_eigen_solve(double g, FockStateSet set);
double* fock_eigenset_get_vector(const EigenSolutionSet set, uint n);
double fock_eigenset_get_value(const EigenSolutionSet set, uint n);
void fock_eigen_free(EigenSolutionSet* es);

/* === Printing and saving to file === */
/* =================================== */

void fock_single_particle_fprint(FILE* f, uint state_nr, StateType type);
static inline void fock_single_particle_print(uint state_nr, StateType type) {
	fock_single_particle_fprint(stdout, state_nr, type);
}

void fock_fprint(FILE* f, const FockState* fock);
static inline void fock_print(const FockState* fock) {
	fock_fprint(stdout, fock);
}

void fock_stateset_fprint(FILE* f, const FockStateSet set);
static inline void fock_stateset_print(const FockStateSet set) {
	fock_stateset_fprint(stdout, set);
}

void fock_stateset_fprint_weighted(FILE* f, const FockStateSet set, const double* coeff, uint n);
static inline void fock_stateset_print_weighted(const FockStateSet set, const double* coeff, uint n) {
	fock_stateset_fprint_weighted(stdout, set, coeff, n);
}

void fock_coeffstateset_fprint(FILE* f, const FockCoeffStateSet set, uint n);
static inline void fock_coeffstateset_print(const FockCoeffStateSet set, uint n) {
	fock_coeffstateset_fprint(stdout, set, n);
}

void fock_stateset_tofile(const char* filename, const FockStateSet set);
FockStateSet fock_stateset_alloc_fromfile(const char* filename);

void fock_stateset_image_to_file(const char* fname_phase, const char* fname_density, FockStateSet set, const double* coeff, double* x, double *y, uint n);

/* === Some additional utility === */
/* =============================== */

#define QF_ARR_LEN 512
#define QF_DR (1.0/(QF_ARR_LEN-1))

typedef struct {
        double arr[QF_ARR_LEN];
} QuickFunction;

double qf_get_value(const QuickFunction* qf, double r);
double qf_integrate(const QuickFunction* qf);
const QuickFunction* bessel_get_qf(uint nr);

void util_bin_tofile(const void* src, uint nbytes, const char* filename);
void* util_bin_alloc_fromfile(uint nbytes, const char* filename);

typedef struct {
        double* data;
        uint nx;
        uint ny;
} UtilMatrix;

void util_numpy_save(const char* filename, double* data, uint nrow, uint ncol);
UtilMatrix util_numpy_load_alloc(const char* filename);
void util_numpy_matrix_free(UtilMatrix* matrix);

#endif /* _FOCK_H_ */
