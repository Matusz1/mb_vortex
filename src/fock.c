#include "fock.h"
#include "fock_internal.h"
#include "util.h"

#include <string.h>
#include <tgmath.h>
#include <time.h>

#include <omp.h>
#include <lapacke.h>

#include <gsl/gsl_sort_uint.h>
#include <gsl/gsl_heapsort.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>

#define FOCK_CONFIG_SET_WEIGHT_CUTOFF 1.0e-4

/* === Lattice === */

Lattice2D fock_lattice_alloc(Vec2d ll, Vec2d ur, uint n_rows, uint n_cols) {
        util_error(n_rows == 0 || n_cols == 0, "Wrong %s arguments.\n", __func__);
        Lattice2D ret = {
                .x = util_malloc(sizeof(*ret.x) * n_rows * n_cols),
                .y = util_malloc(sizeof(*ret.y) * n_rows * n_cols),
                .size = n_rows * n_cols
        };
        uint ij = 0;
        const double lx = ur.x - ll.x;
        const double ly = ur.y - ll.y;
        const double dx = n_cols == 1 ? 0.0 : lx / (n_cols-1);
        const double dy = n_rows == 1 ? 0.0 : ly / (n_rows-1);
        for (uint i = 0; i < n_rows; ++i) {
                for (uint j = 0; j < n_cols; ++j, ++ij) {
                        ret.x[ij] = ll.x + j * dx;
                        ret.y[ij] = ur.y - i * dy;
                }
        }
        return ret;
}

Lattice2D fock_lattice_alloc_helix(uint n_points, double r) {
        const double dphi = 1.0;
        const double alpha = r / sqrt(n_points*dphi);
        Lattice2D ret = {
                .size = n_points,
                .x = util_malloc(sizeof(*ret.x)*n_points),
                .y = util_malloc(sizeof(*ret.y)*n_points),
        };

        for (uint i = 0; i < n_points; ++i) {
                double phi = i*dphi;
                double r = alpha * sqrt(phi);
                ret.x[i] = r * cos(phi);
                ret.y[i] = r * sin(phi);
        }
        return ret;
}

void fock_lattice_free(Lattice2D* lattice) {
        free(lattice->x);
        free(lattice->y);
        lattice->size = 0;
}

/* === Fock states === */

FockState fock_init(uint size, StateType type) {
	util_error(size > FOCK_CONFIG_STATE_CAPACITY, "fock_init, size too big\n");
        FockState fock = {
                .size = size,
                .type = type,
                .states = { 0 },
        };
	return fock;
}

FockState fock_init_direct(uint size, StateType type, const uint* nrs) {
	FockState fock = fock_init(size, type);
	memcpy(fock.states, nrs, sizeof(*nrs)*size);
	gsl_sort_uint(fock.states, 1, size);
	return fock;
}

FockState fock_init_quick_bessel(uint size, const int* m, const uint* rad_exc) {
        uint nrs[size];
        for (uint i = 0; i < size; ++i)
                nrs[i] = fock_bessel_find_index(fock_bessel_state_direct(m[i], rad_exc[i]));
        return fock_init_direct(size, STATE_BESSEL, nrs);
}

bool fock_equal(const FockState* s1, const FockState* s2) {
	if (s1->type != s2->type)
		return false;
	if (s1->size != s2->size)
		return false;
	if (memcmp(s1->states, s2->states, s1->size*sizeof(*s1->states)) != 0)
		return false;
	return true;
}

static inline bool __is_sorted(uint* arr, uint n) {
	for (uint i = 1; i < n; ++i)
		if (arr[i] < arr[i-1])
			return false;
	return true;
}

double fock_single_energy(uint state_nr, StateType type) {
	switch (type) {
	case STATE_PLANEWAVE:
		return __planewave_energy(state_nr);
        case STATE_BESSEL:
                return __bessel_energy(state_nr);
	default:
		util_error(true, "unknown / unsupported type in function %s\n", __func__);
	}
}

double fock_energy(const FockState* fock) {
	double E = 0.0;
	for (uint i = 0; i < fock->size; ++i)
		E += fock_single_energy(fock->states[i], fock->type);
	return E;
}

/* Calculating the values of wave functions */
/* ======================================== */

static void __dec_to_binarr(int* res, uint n, int dim) {
	int pos = dim - 1;
	for (uint i = 0; i != dim+1; ++i)
		res[i] = 0;

	while (n != 0) {
		res[pos] = n & 1;
		res[dim] += res[pos];
                n = n >> 1;
		//n = n / 2;
		--pos;
	}
}

#define MINUS_ONE_POW(n) (((n)%2)==0 ? 1 : -1)

double complex permanent_z(complex double* A, int n) {
	static int chi[FOCK_CONFIG_STATE_CAPACITY+1];
        uint C = (1 << n);
	//const double C = pow(2, n); 
	complex double sum = 0;
	complex double rowsumprod, rowsum;

	for (uint k = 1; k < C; ++k) {
		rowsumprod = 1;
		__dec_to_binarr(chi, k, n);

		for (int m = 0; m < n; ++m) {
			rowsum = 0;
			for (int p = 0; p < n; ++p)
				rowsum += chi[p] * A[m * n + p];
			rowsumprod *= rowsum;    
		}        

		sum += MINUS_ONE_POW(n-chi[n]) * rowsumprod;
		//sum += pow(-1.0, n - chi[n]) * rowsumprod;
	}    

	return sum;
}

/* === One particle wavefunctions === */

complex double fock_single_compute(uint state_nr, StateType type, double x, double y) {
	switch(type) {
	case STATE_PLANEWAVE:
		return __planewave_value(state_nr, x, y);
        case STATE_BESSEL:
                return __bessel_value_xy(state_nr, x, y);
	default:
		util_error(true, "unknown / unsupported type in function %s\n", __func__);
	}
}

Vec2z __single_particle_compute_grad(uint state_nr, StateType type, double x, double y) {
        switch (type) {
        case STATE_BESSEL:
                return __fock_bessel_grad(state_nr, x, y);
        case STATE_PLANEWAVE:
                return __fock_planewave_grad(state_nr, x, y);
        default:
                util_error(true, "Unsupported state.\n");
                return (Vec2z) { 0 };
        }
}

complex double fock_compute(FockState *fock , double *x, double *y) {
	complex double A[fock->size*fock->size];
	uint ij = 0;
	for (uint i = 0; i < fock->size; ++i)
		for (uint j = 0; j < fock->size; ++j, ++ij)
			A[ij] = fock_single_compute(fock->states[j], fock->type, x[i], y[i]);
        double v = 1.0;
        for (uint i = 1; i <= fock->size; ++i)
                v *= i;
	return permanent_z(A, fock->size) / sqrt(v);
}

 Vec2z fock_compute_grad(FockState* fock, double* x, double* y, uint n) {
	complex double Ax[fock->size*fock->size];
	complex double Ay[fock->size*fock->size];
	uint ij = 0;
	for (uint i = 0; i < fock->size; ++i) {
		for (uint j = 0; j < fock->size; ++j, ++ij) {
                        if (i == n) {
                                Vec2z grad = __single_particle_compute_grad(fock->states[j], fock->type, x[i], y[i]);
                                Ax[ij] = grad.x;
                                Ay[ij] = grad.y;
                        }
                        else {
			        Ax[ij] = fock_single_compute(fock->states[j], fock->type, x[i], y[i]);
                                Ay[ij] = Ax[ij];
                        }
                }
        }
        double v = 1.0;
        for (uint i = 1; i <= fock->size; ++i)
                v *= i;
        Vec2z ret = {
                .x = permanent_z(Ax, fock->size) / sqrt(v),
                .y = permanent_z(Ay, fock->size) / sqrt(v),
        };
	return ret;
}

double fock_density(FockState *fock , double *x, double *y) {
	double v = fabs(fock_compute(fock, x, y));
	return v*v;
}

double fock_phase(FockState *fock , double *x, double *y) {
	return carg(fock_compute(fock, x, y));
}

/* === Way of handling sets of fock states === */
/* =========================================== */

void __fock_stateset_alloc_recursive(FockStateSet* set, FockState* running, uint n_states, uint n_particles, double e_cutoff) {
        if (n_particles == 0) {
                fock_stateset_append(set, running);
                return;
        }
        for (uint i = 0; i <= running->states[n_particles]; ++i) {
                running->states[n_particles-1] = i;
                if (fock_energy(running) <= e_cutoff)
                        __fock_stateset_alloc_recursive(set, running, n_states, n_particles-1, e_cutoff);
        }
        running->states[n_particles-1] = 0; /* Reset to the least energy state */
}

FockStateSet fock_stateset_alloc(StateType type, uint n_states, uint n_particles, double e_cutoff) {
        FockStateSet set = { 0 };
        FockState tmp = {
                .size = n_particles,
                .states = { 0 },
                .type = type
        };
        for (uint i = 0; i < n_states; ++i) {
                tmp.states[n_particles-1] = i;
                if (fock_energy(&tmp) <= e_cutoff)
                        __fock_stateset_alloc_recursive(&set, &tmp, n_states, n_particles-1, e_cutoff);
        }
        return set;
}

void __fock_stateset_remove_state(FockStateSet* set, uint idx) {
	util_error(idx >= set->size, "Trying to remove nonexistent element.\n");
	memcpy(&set->states[idx], &set->states[set->size-1], sizeof(*set->states));
	--set->size;
	set->states = realloc(set->states, sizeof(*set->states)*set->size);
}

void __fock_stateset_remove_duplicates(FockStateSet* set) {
	for (uint i = 0; i < set->size; ++i)
		for (uint j = i+1; j < set->size;)
			if (fock_equal(&set->states[i], &set->states[j]))
				__fock_stateset_remove_state(set, j);
			else
				++j;
}

FockStateSet fock_stateset_alloc_copy(FockStateSet src) {
        FockStateSet ret = {
                .states = util_malloc(sizeof(*src.states)*src.size),
                .size = src.size
        };
        memcpy(ret.states, src.states, sizeof(*src.states)*src.size);
        return ret;
}

void fock_stateset_append(FockStateSet* set, const FockState* fock) {
	++set->size;
	set->states = realloc(set->states, sizeof(*set->states) * set->size);
	util_error(set->states == NULL, "Failed to realloc.\n");
	memcpy(set->states + (set->size - 1), fock, sizeof(*fock));
}

void fock_stateset_free(FockStateSet* set) {
	free(set->states);
	set->size = 0;
}

double* fock_stateset_perturbation_function(FockStateSet left, FockState* state) {
        double* coeff = util_malloc(sizeof(*coeff)*left.size);
        for (uint i = 0; i < left.size; ++i) {
                if (fock_equal(&left.states[i], state)) {
                        coeff[i] = 1.0;
                        continue;
                }
                const double v = fock_matrix_element_delta_potential(&left.states[i], state);
                const double dE = fock_energy(state) - fock_energy(&left.states[i]);
                util_error(v != 0.0 && dE==0, "dE = 0 in perturbation theory\n");
                coeff[i] = v / dE;
        }
        return coeff;
}

double fock_stateset_perturbation_energy_one(double g, FockState *state) {
        return g * fock_matrix_element_delta_potential(state, state);
}

double fock_stateset_perturbation_energy_two(double g, FockStateSet left, FockState *state) {
        double ret = 0.0;
        for (uint i = 0; i < left.size; ++i) {
                if (fock_equal(&left.states[i], state))
                        continue;
                const double v = g * fock_matrix_element_delta_potential(&left.states[i], state);
                double dE = fock_energy(state) - fock_energy(&left.states[i]);
                util_error(v != 0.0 && dE==0, "dE = 0 in perturbation theory\n");
                ret += v*v/dE;
        }
        return ret;
}

double fock_stateset_perturbation_energy_three(double g, FockStateSet set, FockState *state) {
        double ret = 0.0;
        for (uint k = 0; k < set.size; ++k) {
                if (fock_equal(&set.states[k], state))
                        continue;
                double v3 = g * fock_matrix_element_delta_potential(&set.states[k], state);
                if (v3 == 0.0)
                        continue;
                double dE3 = fock_energy(&set.states[k]) - fock_energy(state);

                double sum = 0.0;
                for (uint m = 0; m < set.size; ++m) {
                        if (fock_equal(&set.states[m], state))
                                continue;
                        double v1 = g * fock_matrix_element_delta_potential(state, &set.states[m]);
                        if (v1 == 0.0)
                                continue;
                        double v2 = g * fock_matrix_element_delta_potential(&set.states[m], &set.states[k]);
                        if (v2 == 0.0)
                                continue;
                        double dE1 = fock_energy(&set.states[m]) - fock_energy(state);
                        util_error(dE1 == 0.0 || dE3 == 0.0, "dE = 0 in perturbation theory\n");
                        sum += v1*v2*v3 / (dE1*dE3);
                }
                ret += sum;
        }
        double sum = 0.0;
        for (uint m = 0; m < set.size; ++m) {
                if (fock_equal(&set.states[m], state))
                        continue;
                double v = g * fock_matrix_element_delta_potential(state, &set.states[m]);
                if (v == 0.0)
                        continue;
                double dE = fock_energy(&set.states[m]) - fock_energy(state);
                util_error(dE == 0.0, "dE = 0 in perturbation theory\n");
                sum += v*v / (dE*dE);
        }
        sum *= g * fock_matrix_element_delta_potential(state, state);
        return ret - sum;
}

uint fock_stateset_find(FockStateSet set, const FockState* state) {
	uint i = 0;
	for (; i < set.size; ++i)
		if (fock_equal(state, &set.states[i]))
			break;
	return i;
}

uint fock_operator_count(const FockState* fock, uint nr) {
	uint count = 0;
	for (uint i = 0; i < fock->size; ++i)
		if (fock->states[i] == nr)
			++count;
	return count;
}

double fock_operator_create(FockState* fock, uint nr) {
	util_error(fock->size == FOCK_CONFIG_STATE_CAPACITY, "Cannot add more particles to fock state.\n");
	uint idx = fock->size;
        for (; idx > 0; --idx)
                if (nr < fock->states[idx-1])
                        fock->states[idx] = fock->states[idx-1];
                else
                        break;
	fock->states[idx] = nr;
	++fock->size;
	return sqrt(fock_operator_count(fock, nr));
}

double fock_operator_annihilate(FockState* fock, uint nr) {
	uint idx = 0;
	uint count = 0;
	for (uint i = 0; i < fock->size; ++i) {
		if (nr == fock->states[i]) {
			++count;
			idx = i;
		}
	}
	/* Special case when we cannot annichilate any particles,
	 * state just becomes 0 */
	if (count == 0) {
		fock->size = 0;
		return 0.0;
	}
	memmove(fock->states+idx, fock->states+idx+1, sizeof(*fock->states)*(fock->size-idx-1));
	--fock->size;
	return sqrt((double)count);
}

static void __coeffstateset_append_unchecked(FockCoeffStateSet cset[static 1], FockState fstate[static 1], uint capacity[static 1], complex double coeff) {
        const uint idx = fock_stateset_find(cset->set, fstate);
        /* There is fockstate in set, we have to add new one */
        if (idx == cset->set.size) {
                if (*capacity == cset->set.size) {
                        const uint new_cap = (*capacity) * 2;
                        const uint s1 = sizeof(*cset->coeff)*new_cap;
                        const uint s2 = sizeof(*cset->set.states)*new_cap;
                        cset->coeff = realloc(cset->coeff, s1);
                        cset->set.states = realloc(cset->set.states, s2);
                        *capacity = new_cap;
                }
                cset->set.states[cset->set.size] = *fstate;
                cset->coeff[cset->set.size] = coeff;
                ++cset->set.size;
        } else {
                cset->coeff[idx] += coeff;
        }
}

FockCoeffStateSet fock_operator_psi(FockCoeffStateSet cset, double x, double y) {
        uint max = 0;
        for (uint i = 0; i < cset.set.size; ++i)
                if (max < cset.set.states[i].states[cset.set.states[i].size-1])
                        max = cset.set.states[i].states[cset.set.states[i].size-1];
        ++max;

        uint capacity = 16;
        FockCoeffStateSet ret = {
                .coeff = util_malloc(sizeof(*ret.coeff)*capacity),
                .set = {
                        .size = 0,
                        .states = util_malloc(sizeof(*ret.set.states)*capacity),
                },
        };


        FockState spstate;
        /* For every annichilation operator */
        for (uint i = 0; i < max; ++i) {
                const complex double vpsi = fock_single_compute(i, cset.set.states[0].type, x, y);
                /* Apply the operator to every element of set */
                for (uint j = 0; j < cset.set.size; ++j) {
                        complex double v  = cset.coeff[j];
                        if (FOCK_CONFIG_SET_WEIGHT_CUTOFF < cabs(v)) {
                                spstate = cset.set.states[j];
                                v *= fock_operator_annihilate(&spstate, i) * vpsi;
                                if (v != 0.0)
                                        __coeffstateset_append_unchecked(&ret, &spstate, &capacity, v);
                        }
                }
        }

        ret.coeff = util_realloc(ret.coeff, sizeof(*ret.coeff)*ret.set.size);
        double norm = 0.0;
        for (uint i = 0; i < ret.set.size; ++i) {
                double v = cabs(ret.coeff[i]);
                norm += v*v;
        }
        cblas_zdscal(ret.set.size, 1/sqrt(norm), ret.coeff, 1);
        ret.set.states = util_realloc(ret.set.states, sizeof(*ret.set.states)*ret.set.size);
        return ret;
}

FockCoeffStateSet fock_coeffstateset_alloc_dcoeff(FockStateSet set, double* coeff, double cutoff) {
        FockCoeffStateSet ret = {
                .coeff = util_malloc(sizeof(*ret.coeff)*set.size),
                .set = {
                        .size = 0,
                        .states = util_malloc(sizeof(*ret.set.states)*set.size),
                },
        };
        uint j = 0;
        for (uint i = 0; i < set.size; ++i) {
                if (fabs(coeff[i])*fabs(coeff[i]) > cutoff) {
                        ret.set.states[j] = set.states[i];
                        ret.coeff[j] = coeff[i];
                        ++ret.set.size;
                        ++j;
                }
        }
        ret.coeff = util_realloc(ret.coeff, sizeof(*ret.coeff)*ret.set.size);
        ret.set.states = util_realloc(ret.set.states, sizeof(*ret.set.states)*ret.set.size);
        return ret;
}

void fock_coeffstateset_free(FockCoeffStateSet* set) {
        fock_stateset_free(&set->set);
        free(set->coeff);
}

complex double fock_coeffstateset_compute(FockCoeffStateSet cset, double* x, double* y) {
	complex double ret = 0.0;
        //#pragma omp parallel for
	for (uint i = 0; i < cset.set.size; ++i)
		if (FOCK_CONFIG_SET_WEIGHT_CUTOFF < cabs(cset.coeff[i]))
			ret += cset.coeff[i] * fock_compute(&cset.set.states[i], x, y);
	return ret;
}

double fock_matrix_element_delta_potential(const FockState* left, const FockState* right) {
	util_error(left->type != right->type, "Different types of states while calculating matrix elements\n");
	switch (left->type) {
	case STATE_PLANEWAVE:
		return __fock_planewave_matrix_element_delta_potential(left, right);
        case STATE_BESSEL:
                return __bessel_matrix_element_delta_potential(left, right);
	default:
		util_error(true, "unknown / unsupported type in function %s\n", __func__);
	}
}

double fock_matrix_element_kinetic(const FockState* left, const FockState* right) {
	util_error(left->type != right->type, "Different types of states while calculating matrix elements\n");
	if (!fock_equal(left, right))
		return 0.0;
	return fock_energy(left);
}

complex double fock_stateset_compute_weighted(FockStateSet set, double* x, double* y, const double* coeff) {
	complex double ret = 0.0;
        //#pragma omp parallel for
	for (uint i = 0; i < set.size; ++i)
		if (FOCK_CONFIG_SET_WEIGHT_CUTOFF < fabs(coeff[i]))
			ret += coeff[i] * fock_compute(&set.states[i], x, y);
	return ret;
}

Vec2z fock_stateset_compute_weighted_grad(FockStateSet set, double* x, double* y, const double* coeff, uint n) {
        Vec2z ret = { 0 };
        //#pragma omp parallel for reduction(+:tot_x,tot_y)
	for (uint i = 0; i < set.size; ++i) {
		if (FOCK_CONFIG_SET_WEIGHT_CUTOFF < fabs(coeff[i])) {
                        Vec2z grad = fock_compute_grad(&set.states[i], x, y, n);
			ret.x += coeff[i] * grad.x;
			ret.y += coeff[i] * grad.y;
                }
        }
	return ret;
}

double fock_stateset_compute_density_weighted(FockStateSet set, double* x, double* y, const double* coeff) {
	const double v = fabs(fock_stateset_compute_weighted(set, x, y, coeff));
	return v*v;
}

double fock_stateset_compute_phase_weighted(FockStateSet set, double* x, double* y, const double* coeff) {
	return carg(fock_stateset_compute_weighted(set, x, y, coeff));
}

void fock_stateset_bohm_next(FockStateSet set, double* x, double* y, const double* coeff, double dt) {
        const uint n_part = set.states[0].size;
        double vx[n_part];
        double vy[n_part];
        const double complex psi = fock_stateset_compute_weighted(set, x, y, coeff);
        for (uint i = 0; i < n_part; ++i) {
                Vec2z grad = fock_stateset_compute_weighted_grad(set, x, y, coeff, i);
                vx[i] = cimag(grad.x/psi);
                vy[i] = cimag(grad.y/psi);
        }
        for (uint i = 0; i < n_part; ++i) {
                x[i] += vx[i] * dt;
                y[i] += vy[i] * dt;
        }
}

static void __transpose(double* a, uint size) {
	for (uint i = 0; i < size; ++i) {
		for (uint j = i; j < size; ++j) {
			double tmp = a[i*size + j];
			a[i*size + j] = a[j*size + i];
			a[j*size + i] = tmp;
		}
	}
}

void fock_matrix_fill(FockStateSet set, double* T, double* V) {
        //#pragma omp parallel for collapse(2)
        for (uint i = 0; i < set.size; ++i) {
                for (uint j = i; j < set.size; ++j) {
                        double t = fock_matrix_element_kinetic(&set.states[i], &set.states[j]);
                        double v = fock_matrix_element_delta_potential(&set.states[i], &set.states[j]);
                        V[i*set.size + j] = v;
                        V[j*set.size + i] = v;
                        T[i*set.size + j] = t;
                        T[j*set.size + i] = t;
                }
        }
}

EigenSolutionSet fock_eigen_vt_solve(double g, uint n, double* T, double* V) {
        EigenSolutionSet set = {
                .size = n,
                .vectors = util_malloc(sizeof(*set.vectors)*n*n),
                .values = util_malloc(sizeof(*set.values)*n)
        };
        for (uint i = 0; i < n; ++i) {
                for (uint j = i; j < n; ++j) {
                        set.vectors[i*n+j] = T[i*n+j] + g*V[i*n+j];
                }
        }
	LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', set.size, set.vectors, set.size, set.values);
	__transpose(set.vectors, set.size);
	return set;
}

EigenSolutionSet fock_eigen_omega_vt_solve(double g, double omega, FockStateSet fset, double* T, double* V) {
        EigenSolutionSet set = {
                .size = fset.size,
                .vectors = util_malloc(sizeof(*set.vectors)*fset.size*fset.size),
                .values = util_malloc(sizeof(*set.values)*fset.size)
        };
        for (uint i = 0; i < fset.size; ++i) {
                for (uint j = i; j < fset.size; ++j) {
                        set.vectors[i*fset.size+j] = T[i*fset.size+j] + g*V[i*fset.size+j];
                }
        }
        for (uint i = 0; i < set.size; ++i) {
                set.vectors[i*set.size+i] -= omega * fock_bessel_angular_momentum(&fset.states[i]);
        }
	LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', set.size, set.vectors, set.size, set.values);
	__transpose(set.vectors, set.size);
	return set;
}

EigenSolutionSet fock_eigen_solve(double g, FockStateSet set) {
	EigenSolutionSet esols = {
		.size = set.size,
		.vectors = util_malloc(sizeof(*esols.vectors)*set.size*set.size),
		.values = util_malloc(sizeof(*esols.values)*set.size)
	};
        //#pragma omp parallel for collapse(2)
	for (uint i = 0; i < set.size; ++i) {
		for (uint j = i; j < set.size; ++j) {
			const double t = fock_matrix_element_kinetic(&set.states[i], &set.states[j]);
			const double v = fock_matrix_element_delta_potential(&set.states[i], &set.states[j]);
			esols.vectors[i*set.size + j] = t + g*v;
		}
	}
	LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', set.size, esols.vectors, set.size, esols.values);
	__transpose(esols.vectors, esols.size);
	return esols;
}

double* fock_eigenset_get_vector(const EigenSolutionSet set, uint n) {
        util_error(set.size <= n, "Too big of an index.\n");
        return set.vectors + n*set.size;
}

double fock_eigenset_get_value(const EigenSolutionSet set, uint n) {
        util_error(set.size <= n, "Too big of an index.\n");
        return set.values[n];
}

void fock_eigen_free(EigenSolutionSet* es) {
	free(es->vectors);
	free(es->values);
}

void fock_single_particle_fprint(FILE* f, uint state_nr, StateType type) {
	switch (type) {
	case STATE_PLANEWAVE:
		__fock_planewave_fprint(f, state_nr);
		break;
        case STATE_BESSEL:
                __fock_bessel_fprint(f, state_nr);
                break;
	default:
		util_error(true, "unknown / unsupported type in function %s\n", __func__);
		break;
	}
}

void fock_fprint(FILE* f, const FockState* fock) {
	for (uint i = 0; i < fock->size;) {
		fock_single_particle_fprint(f, fock->states[i], fock->type);
		if (++i != fock->size)
			putc(' ', f);
	}
}

void fock_stateset_fprint(FILE* f, const FockStateSet set) {
	for (uint i = 0; i < set.size; ++i) {
		fprintf(f, "%6u: ", i);
		fock_fprint(f, &set.states[i]);
		putc('\n', f);
	}
}

static int __fock_square_compare(const void* l, const void* r) {
	double al = fabs(*(double*)l);
	double ar = fabs(*(double*)r);
	if (al > ar)
		return -1;
	else if (al == ar)
		return 0;
	else
		return 1;
}

static int __fock_complex_square_compare(const void* l, const void* r) {
	double al = cabs(*(complex double*)l);
	double ar = cabs(*(complex double*)r);
	if (al > ar)
		return -1;
	else if (al == ar)
		return 0;
	else
		return 1;
}

void fock_stateset_fprint_weighted(FILE* f, const FockStateSet set, const double* coeff, uint n) {
	size_t idxs[set.size];
	gsl_heapsort_index(idxs, coeff, set.size, sizeof(*coeff), __fock_square_compare);
	fprintf(f, " |ci|^2         state\n");
	for (uint i = 0; i < MIN(set.size, n); ++i) {
		fprintf(f, "%lf : ", coeff[idxs[i]]*coeff[idxs[i]]);
		fock_fprint(f, &set.states[idxs[i]]);
		putc('\n', f);
	}
}

void fock_coeffstateset_fprint(FILE* f, const FockCoeffStateSet cset, uint n) {
	size_t idxs[cset.set.size];
	gsl_heapsort_index(idxs, cset.coeff, cset.set.size, sizeof(*cset.coeff), __fock_complex_square_compare);
	fprintf(f, " |ci|^2         state\n");
	for (uint i = 0; i < MIN(cset.set.size, n); ++i) {
		fprintf(f, "%lf : ", cabs(cset.coeff[idxs[i]])*cabs(cset.coeff[idxs[i]]));
		fock_fprint(f, &cset.set.states[idxs[i]]);
		putc('\n', f);
	}
}

#define __FOCK_STATESET_CONTROL_NUM 1983456

void fock_stateset_tofile(const char* filename, const FockStateSet set) {
	const uint32_t ctr_num = __FOCK_STATESET_CONTROL_NUM;
	FILE* f = fopen(filename, "wb");
	util_error(f == NULL, "Failed to open file, while writing stateset.\n");
	util_fwrite(&ctr_num, sizeof(ctr_num), 1, f);
	util_fwrite(&set.size, sizeof(set.size), 1, f);
	util_fwrite(set.states, sizeof(*set.states), set.size, f);
	fclose(f);
}

FockStateSet fock_stateset_alloc_fromfile(const char* filename) {
	FILE* f = fopen(filename, "rb");
	util_error(f == NULL, "Failed to open file %s, while writing stateset.\n", filename);
	uint32_t ctr_num;
	util_fread(&ctr_num, sizeof(ctr_num), 1, f);
	util_error(ctr_num != __FOCK_STATESET_CONTROL_NUM, "Control numbers didn't mach wchile reading stateset from %s\n", filename);
	FockStateSet set;
	util_fread(&set.size, sizeof(set.size), 1, f);
	set.states = util_malloc(sizeof(*set.states) * set.size);
	util_fread(set.states, sizeof(*set.states), set.size, f);
	fclose(f);
	return set;
}

static inline double __period(double x) {
	double tmp;
	x = modf(x, &tmp);
	return (x < 0.0) ? (x + 1.0) : x;
}

double* __fock_stateset_planewave_metropolis_alloc(FockStateSet set, const double* coeff, uint n_draws, uint n_burnout, unsigned long seed) {
        const uint n_particles = set.states[0].size;
        double* ret = util_malloc(sizeof(*ret)*n_draws*2*n_particles);
        double x[n_particles];
        double y[n_particles];
        double x_next[n_particles];
        double y_next[n_particles];

        const double sigma = 0.1;

	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, seed);

	for (uint i = 0; i < n_particles; ++i) {
                x[i] = gsl_rng_uniform(rng);
                y[i] = gsl_rng_uniform(rng);
	}

        uint ret_idx = 0;
	for (uint i = 0; i < n_burnout + n_draws; ++i) {
                printf("%5u / %5u\n", i+1, (n_burnout+n_draws));
		for (uint j = 0; j != n_particles; ++j) {
                        x_next[j] = x[j] + gsl_ran_gaussian(rng, sigma);
                        y_next[j] = y[j] + gsl_ran_gaussian(rng, sigma);
		}

		const double p_old = fock_stateset_compute_density_weighted(set, x, y, coeff);
		const double p_new = fock_stateset_compute_density_weighted(set, x_next, y_next, coeff);
		const double alpha = p_new / p_old;
		const double r = gsl_rng_uniform(rng);

		if (r <= alpha) {
			for (uint j = 0; j < n_particles; ++j) {
				x[j] = x_next[j];
				y[j] = y_next[j];
			}
		}

		if (n_burnout <= i) {
			for (uint j = 0; j < n_particles; ++j) {
                                ret[ret_idx] = __period(x[j]);
                                ++ret_idx;
			}
			for (uint j = 0; j < n_particles; ++j) {
                                ret[ret_idx] = __period(y[j]);
                                ++ret_idx;
			}
		}

	}

	gsl_rng_free(rng);
        return ret;
}

double fock_stateset_aiaj(FockStateSet set, const double* coeff, uint i, uint j) {
        FockStateSet set_copy = fock_stateset_alloc_copy(set);
        double op_val[set_copy.size];
        for (uint n = 0; n < set_copy.size; ++n) {
                op_val[n]  = fock_operator_annihilate(set_copy.states+n, j);
                op_val[n] *= fock_operator_create(set_copy.states+n, i);
        }
        double ret = 0.0;
        for (uint n = 0; n < set.size; ++n) {
                if (op_val[n] == 0.0)
                        continue;

                for (uint m = 0; m < set.size; ++m) {
                        if (fock_equal(&set_copy.states[n], &set.states[m])) {
                                ret += coeff[n]*coeff[m]*op_val[n];
                                break;
                        }
                }
        }
        fock_stateset_free(&set_copy);
        return ret;
}

complex double fock_coeffstateset_aiaj(FockCoeffStateSet cset, uint i, uint j) {
        FockStateSet set_copy = fock_stateset_alloc_copy(cset.set);
        double op_val[set_copy.size];
        for (uint n = 0; n < set_copy.size; ++n) {
                op_val[n]  = fock_operator_annihilate(set_copy.states+n, j);
                op_val[n] *= fock_operator_create(set_copy.states+n, i);
        }
        complex double ret = 0.0;
        for (uint n = 0; n < cset.set.size; ++n) {
                if (op_val[n] == 0.0)
                        continue;

                for (uint m = 0; m < cset.set.size; ++m) {
                        if (fock_equal(&set_copy.states[n], &cset.set.states[m])) {
                                ret += cset.coeff[n]*conj(cset.coeff[m])*op_val[n];
                                break;
                        }
                }
        }
        fock_stateset_free(&set_copy);
        return ret;
}


G1State fock_g1_alloc(FockStateSet set, const double* coeff) {
        uint max_ij = 0;
        for (uint i = 0; i < set.size; ++i)
                if (max_ij < set.states[i].states[set.states[i].size-1])
                        max_ij = set.states[i].states[set.states[i].size-1];
        ++max_ij;
        G1State g1 = {
                .size = max_ij*max_ij,
                .type = set.states[0].type,
                .is = util_malloc(sizeof(*g1.is)*max_ij*max_ij),
                .js = util_malloc(sizeof(*g1.js)*max_ij*max_ij),
                .coeff = util_malloc(sizeof(*g1.coeff)*max_ij*max_ij),
        };

        uint ij = 0;
        for (uint i = 0; i < max_ij; ++i) {
                for (uint j = 0; j < max_ij; ++j) {
                        g1.is[ij] = i;
                        g1.js[ij] = j;
                        g1.coeff[ij] = fock_stateset_aiaj(set, coeff, i, j);
                        if (cabs(g1.coeff[ij]) < 1.0e-8)
                                --g1.size;
                        else
                                ++ij;
                }
        }
        g1.coeff = util_realloc(g1.coeff, sizeof(*g1.coeff)*g1.size);
        g1.is = util_realloc(g1.is, sizeof(*g1.is)*g1.size);
        g1.js = util_realloc(g1.js, sizeof(*g1.js)*g1.size);
        return g1;
}

G1State fock_g1_alloc_cset(FockCoeffStateSet cset) {
        uint max_ij = 0;
        for (uint i = 0; i < cset.set.size; ++i)
                if (max_ij < cset.set.states[i].states[cset.set.states[i].size-1])
                        max_ij = cset.set.states[i].states[cset.set.states[i].size-1];
        ++max_ij;
        G1State g1 = {
                .size = max_ij*max_ij,
                .type = cset.set.states[0].type,
                .is = util_malloc(sizeof(*g1.is)*max_ij*max_ij),
                .js = util_malloc(sizeof(*g1.js)*max_ij*max_ij),
                .coeff = util_malloc(sizeof(*g1.coeff)*max_ij*max_ij),
        };

        uint ij = 0;
        for (uint i = 0; i < max_ij; ++i) {
                for (uint j = 0; j < max_ij; ++j) {
                        g1.is[ij] = i;
                        g1.js[ij] = j;
                        g1.coeff[ij] = fock_coeffstateset_aiaj(cset, i, j);
                        if (cabs(g1.coeff[ij]) < 1.0e-8)
                                --g1.size;
                        else
                                ++ij;
                }
        }
        g1.coeff = util_realloc(g1.coeff, sizeof(*g1.coeff)*g1.size);
        g1.is = util_realloc(g1.is, sizeof(*g1.is)*g1.size);
        g1.js = util_realloc(g1.js, sizeof(*g1.js)*g1.size);
        return g1;
}

double* fock_coeffstateset_alloc_random(FockCoeffStateSet cset, uint n) {
        uint n_part = cset.set.states[0].size;
        double* p = util_malloc(sizeof(*p)*n*2*n_part);

        /* Draw first positions */
        fprintf(stderr, "# Allocating first G1 function.\n");
        G1State g6 = fock_g1_alloc_cset(cset);
        fprintf(stderr, "# Drawing first positions\n");
        fock_g1_rand(&g6, FOCK_G1_RAND_RELOAD);
        for (uint i = 0; i < n; ++i) {
                Vec2d v = fock_g1_rand(&g6, FOCK_G1_RAND_KEEP_ARR);
                p[i*(2*n_part)] = v.x;
                p[i*(2*n_part)+n_part] = v.y;
        }
        fock_g1_free(&g6);

        if (n_part == 1)
                return p;

        G1State gfuncs[n_part-1];
        FockCoeffStateSet csets[n_part];
        csets[0] = cset;
        fprintf(stderr, "# Drawing rest of positions.\n");
        for (uint i = 0; i < n; ++i) {
                fprintf(stderr, "[ %4u / %4u ] -> ", i+1, n);
                for (uint j = 0; j < n_part-1; ++j) {
                        printf("%u", j+1);
                        csets[j+1] = fock_operator_psi(csets[j], p[i*(2*n_part)+j], p[i*(2*n_part)+n_part+j]);
                        gfuncs[j] = fock_g1_alloc_cset(csets[j+1]);
                        Vec2d v = fock_g1_rand(gfuncs+j, FOCK_G1_RAND_RELOAD);
                        p[i*(2*n_part)+j+1] = v.x;
                        p[i*(2*n_part)+n_part+j+1] = v.y;
                }
                putchar('\n');


                for (uint j = 0; j <n_part-1; ++j) {
                        fock_coeffstateset_free(csets+j+1);
                        fock_g1_free(gfuncs+j);
                }
        }

        return p;
}

double complex bessel_g1_compute_rphi(G1State* g1, double r, double phi, double rp, double phip) {
        double complex ret = 0.0;
        //#pragma omp parallel for reduction(+:ret)
        for (uint n = 0; n < g1->size; ++n) {
                double complex v = conj(bessel_value(g1->is[n], r, phi));
                v *= bessel_value(g1->js[n], rp, phip);
                ret += g1->coeff[n]*v;
        }
        return ret;
}

double complex fock_g1_compute(G1State* g1, double x, double y, double xp, double yp) {
        double complex ret = 0.0;
        //#pragma omp parallel for reduction(+:ret)
        for (uint n = 0; n < g1->size; ++n) {
                double complex v = conj(fock_single_compute(g1->is[n], g1->type, x, y));
                v *= fock_single_compute(g1->js[n], g1->type, xp, yp);
                ret += g1->coeff[n]*v;
        }
        return ret;
}

static uint __bsearch_double(double key, double* arr, uint n) {
        for (uint i = 0; i < n; ++i)
                if (key <= arr[i])
                        return i;

        return n-1;

        // This binary search is not working,
        // but this unimportant, it would not be quicker
        uint low = 0;
        uint high = n;
        while ((high-low) > 1) {
                uint mid = (high+low)/2;
                if (key < arr[mid])
                        high = mid;
                else
                        low = mid;
        }
        return low;
}

/* Based of CDF (may be unnormalized) raturns values from [0, len-1] */
static inline uint __arr_rand(double* arr, uint len) {
        const double r = arr[len-1] * rand() / RAND_MAX;
        return __bsearch_double(r, arr, len);
}

#define _G1_NR   40
#define _G1_NPHI 40
static double integral_r[_G1_NR];
static double integrals_phi[_G1_NPHI*_G1_NR];

Vec2d fock_g1_rand(G1State* g1, FockG1RandOpt opt) {
        util_error(g1->type != STATE_BESSEL, "Drawing from g1 other than bessel state is not supported.\n");
        const double dr   = 1.0 / _G1_NR;
        const double dphi = 2.0 * M_PI / _G1_NPHI;

        /* Only if proper option is set, we recalculate the integral */
        if (opt == FOCK_G1_RAND_RELOAD) {
                double rsum = 0.0;
                for (uint i = 0; i < _G1_NR; ++i) {
                        const double r = (i+0.5)*dr; 
                        double sum = 0.0;
                        for (uint j = 0; j < _G1_NPHI; ++j) {
                                const double phi = (j+0.5)*dphi; 
                                sum += creal(bessel_g1_compute_rphi(g1, r, phi, r, phi));
                                integrals_phi[i*_G1_NR+j] = sum;
                        }
                        rsum += sum;
                        integral_r[i] = rsum;
                }
        }

        /* Drawing random point */
        uint nr = __arr_rand(integral_r, _G1_NR);
        uint nphi = __arr_rand(&integrals_phi[nr*_G1_NR], _G1_NPHI);
        double r   =   dr * (nr   + ((double)rand()/RAND_MAX));
        double phi = dphi * (nphi + ((double)rand()/RAND_MAX));

        Vec2d ret = {
                .x = r*cos(phi),
                .y = r*sin(phi),
        };

        return ret;
}

#if 0

#define _G1_NX 40
#define _G1_NY 40
double integral[_G1_NX*_G1_NY+1];

Vec2d __old_fock_g1_rand_(G1State* g1, FockG1RandOpt opt) {
        util_error(g1->type != STATE_BESSEL, "Drawing from g1 other than bessel state is not supported.\n");
        const double dx = 2.0 / _G1_NX;
        const double dy = 2.0 / _G1_NY;
        const double dxh = 0.5*dx;
        const double dyh = 0.5*dy;
        const uint ncells = _G1_NX*_G1_NY;
        integral[0] = 0.0;
        uint ij = 0;
        if (opt == FOCK_G1_RAND_RELOAD) {
                for (uint i = 0; i < _G1_NY; ++i) {
                        for (uint j = 0; j < _G1_NX; ++j) {
                                const double x = -1.0 + j*dx + dxh; 
                                const double y = -1.0 + i*dy + dyh; 
                                integral[ij+1] = integral[ij] + creal(fock_g1_compute(g1, x, y, x, y));
                                ++ij;
                        }
                }
        }
        util_numpy_save("lol.dat", integral, _G1_NY, _G1_NX);

        /* Drawing random point */
        double r = integral[ncells] * rand() / RAND_MAX;
        ij = __bsearch_double(r, integral, ncells);
        double diff = integral[ij+1] - integral[ij];
        double ddx = 0.0;
        if (diff > 0.0)
                ddx = dx * (r - integral[ij])/diff;
        double ddy = dy * rand() / RAND_MAX;

        Vec2d ret = {
                .x = -1.0 + dx*(ij % _G1_NX) + ddx,
                .y = -1.0 + dy*((double)(ij / _G1_NX)) + ddy,
        };
        double dot = ret.x*ret.x + ret.y*ret.y;
        if (dot > 1.0) {
                return fock_g1_rand(g1, FOCK_G1_RAND_KEEP_ARR);
        }

        return ret;
}

#endif

void fock_g1_free(G1State* g1) {
        free(g1->is);
        free(g1->js);
        free(g1->coeff);
}

complex double* fock_g1_matrix_alloc(Lattice2D lattice, G1State* g1) {
        complex double* M = util_malloc(sizeof(*M)*lattice.size*lattice.size);

        uint ixy = 0;
        for (uint i = 0; i < lattice.size; ++i)
                for (uint j = 0; j < lattice.size; ++j, ++ixy)
                        M[ixy] = fock_g1_compute(g1, lattice.x[j], lattice.y[j], lattice.x[i], lattice.y[i]);

        double norm_const = 0.0;
        for (uint i = 0; i < lattice.size; ++i)
                norm_const += creal(M[i*(lattice.size+1)]);
        norm_const = 1 / norm_const;
        cblas_zdscal(lattice.size*lattice.size, norm_const, M, 1);
        return M;
}

double* fock_g1_matrix_calc_top_eigenvector(double complex* M, uint n, uint n_iter, double* eval) {
        complex double* vc = util_malloc(sizeof(*vc)*n);
        complex double* vc_copy = util_malloc(sizeof(*vc)*n);
        srand(time(NULL));
        for (uint i = 0; i < n; ++i) {
                double re = (double) rand() / RAND_MAX;
                double im = (double) rand() / RAND_MAX;
                vc[i] = re + I*im;
        }

        const complex double alpha = 1.0;
        const complex double beta = 0.0;
        double invnorm = 1.0;
        for (uint i = 0; i < n_iter; ++i) {
                cblas_zhemv(CblasRowMajor, CblasUpper, n, &alpha, M, n, vc, 1, &beta, vc_copy, 1);
                invnorm = 1 / cblas_dznrm2(n, vc_copy, 1);
                cblas_zdscal(n, invnorm, vc_copy, 1);
                memcpy(vc, vc_copy, sizeof(*vc)*n);
        }
        *eval = 1 / invnorm;

        double* v = util_malloc(sizeof(*v)*n);
        for (uint i = 0; i < n; ++i)
                v[i] = cabs(vc[i]);
        util_free(vc);

        return v;
}

CmplxEigenSolutionSet fock_g1_matrix_solve(double complex* M, uint dim) {
        CmplxEigenSolutionSet ret = {
                .size = dim,
                .vectors = util_malloc(sizeof(*ret.vectors)*dim*dim),
                .values = util_malloc(sizeof(*ret.values)*dim),
        };
        LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', ret.size, M, ret.size, ret.values);
        for (uint i = 0; i < ret.size; ++i)
                for (uint j = 0; j < ret.size; ++j)
                        ret.vectors[i*ret.size+j] = M[j*ret.size+i];

        return ret;
}

void fock_eigen_save(const char* fname, EigenSolutionSet set) {
        FILE* f = util_fopen(fname, "wb");
        util_fwrite(&set.size, sizeof(set.size), 1, f);
        util_fwrite(set.values, sizeof(*set.values), set.size, f);
        util_fwrite(set.vectors, sizeof(*set.vectors), set.size*set.size, f);
        fclose(f);
}

EigenSolutionSet fock_eigen_alloc_file(const char* fname) {
        FILE* f = util_fopen(fname, "rb");
        EigenSolutionSet set;
        util_fread(&set.size, sizeof(set.size), 1, f);
        set.values = util_malloc(sizeof(*set.values)*set.size);
        set.vectors = util_malloc(sizeof(*set.vectors)*set.size*set.size);
        util_fread(set.values, sizeof(*set.values), set.size, f);
        util_fread(set.vectors, sizeof(*set.vectors), set.size*set.size, f);
        fclose(f);
        return set;
}

double* __fock_stateset_bessel_metropolis_alloc(FockStateSet set, const double* coeff, uint n_draws, uint n_burnout, unsigned long seed) {
        const uint n_particles = set.states[0].size;
        double* ret = util_malloc(sizeof(*ret)*n_draws*2*n_particles);
        double x[n_particles];
        double y[n_particles];
        double x_next[n_particles];
        double y_next[n_particles];

        const double sigma = 0.1;

	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, seed);

	for (uint i = 0; i < n_particles; ++i) {
                do {
                        x[i] = gsl_rng_uniform(rng)*2.0-1.0;
                        y[i] = gsl_rng_uniform(rng)*2.0-1.0;
                } while ((x[i]*x[i] + y[i]*y[i]) > 1);
	}

        uint ret_idx = 0;
	for (uint i = 0; i < n_burnout + n_draws; ++i) {
                printf("%u / %u\n", i+1, n_draws+n_burnout);
		for (uint j = 0; j != n_particles; ++j) {
                        do {
                                x_next[j] = x[j] + gsl_ran_gaussian(rng, sigma);
                                y_next[j] = y[j] + gsl_ran_gaussian(rng, sigma);
                        } while ((x_next[j]*x_next[j] + y_next[j]*y_next[j]) > 1);
		}

		const double p_old = fock_stateset_compute_density_weighted(set, x, y, coeff);
		const double p_new = fock_stateset_compute_density_weighted(set, x_next, y_next, coeff);
		const double alpha = p_new / p_old;
		const double r = gsl_rng_uniform(rng);

		if (r <= alpha) {
			for (uint j = 0; j < n_particles; ++j) {
				x[j] = x_next[j];
				y[j] = y_next[j];
			}
		}

		if (n_burnout <= i) {
			for (uint j = 0; j < n_particles; ++j) {
                                ret[ret_idx] = x[j];
                                ++ret_idx;
			}
			for (uint j = 0; j < n_particles; ++j) {
                                ret[ret_idx] = y[j];
                                ++ret_idx;
			}
		}

	}

	gsl_rng_free(rng);
        return ret;

}

double* fock_stateset_metropolis_alloc(FockStateSet set, const double* coeff, uint n_draws, uint n_burnout, unsigned long seed) {
        switch (set.states[0].type) {
        case STATE_PLANEWAVE:
                return __fock_stateset_planewave_metropolis_alloc(set, coeff, n_draws, n_burnout, seed);
        case STATE_BESSEL:
                return __fock_stateset_bessel_metropolis_alloc(set, coeff, n_draws, n_burnout, seed);
        default:
		util_error(true, "unknown / unsupported type in function %s\n", __func__);
        }
}

void util_bin_tofile(const void* src, uint nbytes, const char* filename) {
	FILE* f = util_fopen(filename, "wb");
	util_fwrite(src, 1, nbytes, f);
	fclose(f);
}

void* util_bin_alloc_fromfile(uint nbytes, const char* filename) {
	FILE* f = fopen(filename, "rb");
	util_error(f == NULL, "Failed to open file %s, while reading stateset.\n", filename);
	void* dst = util_malloc(nbytes);
	util_fread(dst, 1, nbytes, f);
	fclose(f);
	return dst;
}

double qf_get_value(const QuickFunction* qf, double r) {
        r = fabs(r);
        if (r > 1.0)
                return 0.0;
        uint l = (uint) (r / QF_DR);
        uint u = l + 1;
        double vl = qf->arr[l];
        double vu = qf->arr[u];
        double dfdr = (vu - vl) / QF_DR;
        double diff = (r - l*QF_DR) * dfdr;
        return vl + diff;
}

double qf_integrate(const QuickFunction* qf) {
        double sum = 0.0;
        for (uint i = 1; i < QF_ARR_LEN-1; ++i) {
                double r = i*QF_DR;
                sum += r*qf->arr[i];
        }
        sum += 0.5*qf->arr[QF_ARR_LEN-1];
        return sum*QF_DR;
}

void __fock_stateset_planewave_image_to_file(const char* fname_phase, const char* fname_density, FockStateSet set, const double* coeff, double* x, double *y, uint n) {
        const double d = 1.0 / n;
        double xp = 0.0;
        double yp = 0.0;
        const double x0 = x[0];
        const double y0 = y[0];
        FILE* f_phase = util_fopen(fname_phase, "w");
        FILE* f_density = util_fopen(fname_density, "w");
        for (uint i = 0; i < n; ++i, yp+=d) {
                xp = 0.0;
                for (uint j = 0; j < n; ++j, xp+=d) {
                        x[0] = xp;
                        y[0] = yp;
                        double complex val = fock_stateset_compute_weighted(set, x, y, coeff);
                        fprintf(f_phase, "%lf ", carg(val));
                        fprintf(f_density, "%lf ", cabs(val)*cabs(val));
                }
                fputc('\n', f_phase);
                fputc('\n', f_density);
        }
        fclose(f_phase);
        fclose(f_density);
        x[0] = x0;
        y[0] = y0;
}

void __fock_stateset_bessel_image_to_file(const char* fname_phase, const char* fname_density, FockStateSet set, const double* coeff, double* x, double *y, uint n) {
        const double d = 2.0 / n;
        double xp = -1.0;
        double yp = -1.0;
        const double x0 = x[0];
        const double y0 = y[0];
        FILE* f_phase = util_fopen(fname_phase, "w");
        FILE* f_density = util_fopen(fname_density, "w");
        for (uint i = 0; i < n; ++i, yp+=d) {
                xp = -1.0;
                for (uint j = 0; j < n; ++j, xp+=d) {
                        x[0] = xp;
                        y[0] = yp;
                        double complex val = 0.0;
                        if ((xp*xp+yp*yp) <= 1.0)
                                val = fock_stateset_compute_weighted(set, x, y, coeff);
                        fprintf(f_phase, "%lf ", carg(val));
                        fprintf(f_density, "%lf ", cabs(val)*cabs(val));
                }
                fputc('\n', f_phase);
                fputc('\n', f_density);
        }
        fclose(f_phase);
        fclose(f_density);
        x[0] = x0;
        y[0] = y0;

}

void fock_stateset_image_to_file(const char* fname_phase, const char* fname_density, FockStateSet set, const double* coeff, double* x, double *y, uint n) {
        switch (set.states[0].type) {
        case STATE_PLANEWAVE:
                __fock_stateset_planewave_image_to_file(fname_phase, fname_density, set, coeff, x, y, n);
                return;
        case STATE_BESSEL:
                __fock_stateset_bessel_image_to_file(fname_phase, fname_density, set, coeff, x, y, n);
                return;
        default:
                util_error(true, "Unknown fock state type in function %s.\n", __func__);
        }
}

#define FOCK_CONFIG_MAX_NUMPY_HEADER_LEN 256

static const char _numpy_magic_str[6] = "\x93NUMPY";

void util_numpy_save(const char* filename, double* data, uint nrow, uint ncol) {
        char header[FOCK_CONFIG_MAX_NUMPY_HEADER_LEN];
        const u8 ver_major = 1;
        const u8 ver_minor = 0;
        sprintf(header, "{'descr': '<f8', 'fortran_order': False, 'shape': (%u, %u), }\n", nrow, ncol);
        u16 header_len = strlen(header);
        u8 spaces_padding_len = 64 - ((header_len+10)%64);
        memset(header+header_len-1, 0x20, spaces_padding_len); /* Padding with spaces */
        header_len += spaces_padding_len;
        header[header_len-1] = '\n';

        FILE* f = util_fopen(filename, "wb");
        util_fwrite(_numpy_magic_str, sizeof(_numpy_magic_str), 1, f);
        util_fwrite(&ver_major, 1, 1, f);
        util_fwrite(&ver_minor, 1, 1, f);
        util_fwrite(&header_len, 2, 1, f);
        util_fwrite(header, header_len, 1, f);
        util_fwrite(data, sizeof(*data), nrow*ncol, f);
        util_fclose(f);
}

UtilMatrix util_numpy_load_alloc(const char* filename) {
        fprintf(stderr, "# WARNING: You are using untested function %s, good luck.\n", __func__);
        char magic_string[6];
        u8 major;
        u8 minor;
        u16 header_len;
        char header[FOCK_CONFIG_MAX_NUMPY_HEADER_LEN];
        FILE* f = util_fopen(filename, "rb");
        UtilMatrix mat;
        //memset(header+header_len+1, 0, FOCK_CONFIG_NUMPY_PADDING_LEN);
        util_fread(magic_string, sizeof(magic_string), 1, f);
        util_error(memcmp(magic_string, _numpy_magic_str, sizeof(_numpy_magic_str)) != 0, "Cannot read numpy file.\n");
        util_fread(&major, sizeof(major), 1, f);
        util_fread(&minor, sizeof(minor), 1, f);
        util_fread(&header_len, sizeof(header_len), 1, f);
        //rgc_fread(header, header_len+FOCK_CONFIG_NUMPY_PADDING_LEN, 1, f);
        util_error(sscanf(header, "(%u, %u)", &mat.nx, &mat.ny) != 2, "Failed to read numpy dimension.\n");
        mat.data = util_malloc(sizeof(*mat.data)*mat.nx*mat.ny);
        util_fread(mat.data, sizeof(*mat.data), mat.nx*mat.ny, f);
        util_fclose(f);
        return mat;
}

void util_numpy_matrix_free(UtilMatrix* matrix) {
        free(matrix->data);
        matrix->nx = 0;
        matrix->ny = 0;
}
