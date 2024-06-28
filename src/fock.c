
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

#include "fock.h"
#include "state.h"
#include "util.h"
#include "rho.h"

#define FOCK_CONFIG_SET_WEIGHT_CUTOFF 1.0e-4

/* === Fock states === */

Fock fock_init(uint size) {
	util_error(size > FOCK_CONFIG_STATE_CAPACITY, "fock_init, size too big\n");
        Fock fock = {
                .size = size,
                .states = { 0 },
        };
	return fock;
}

Fock fock_init_direct(uint size, const uint* nrs) {
	Fock fock = fock_init(size);
	memcpy(fock.states, nrs, sizeof(*nrs)*size);
	gsl_sort_uint(fock.states, 1, size);
	return fock;
}

Fock fock_init_quick(uint size, const int* m, const uint* rad_exc) {
        uint nrs[size];
        for (uint i = 0; i < size; ++i)
                nrs[i] = state_find_index(state_info_direct(m[i], rad_exc[i]));
        return fock_init_direct(size, nrs);
}

void state_basis_free(Basis* basis) {
        util_free(basis->states);
        basis->size = 0;
        basis->capacity = 0;
}

bool state_basis_contains(Basis* basis, uint state) {
        for (uint i = 0; i < basis->size; ++i)
                if (basis->states[i] == state)
                        return true;
        return false;
}

Basis state_basis_alloc_copy(Basis* basis) {
        Basis ret = {
                .states = util_malloc(sizeof(*ret.states)*basis->size),
                .size = basis->size,
                .capacity = basis->size,
        };
        memcpy(ret.states, basis->states, sizeof(*ret.states)*ret.size);
        return ret;
}

void state_basis_add(Basis* basis, uint state) {
        if (basis->size == basis->capacity) {
                basis->capacity = MAX(1, 2*basis->capacity);
                basis->states = util_realloc(basis->states, sizeof(*basis->states)*basis->capacity);
        }
        basis->states[basis->size] = state;
        ++basis->size;
}

bool fock_equal(const Fock* s1, const Fock* s2) {
	if (s1->size != s2->size)
		return false;
	if (memcmp(s1->states, s2->states, s1->size*sizeof(*s1->states)) != 0)
		return false;
	return true;
}

double fock_energy(const Fock* fock) {
	double E = 0.0;
	for (uint i = 0; i < fock->size; ++i)
		E += state_energy(fock->states[i]);
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

complex double fock_compute(Fock *fock , double *x, double *y) {
        const uint s = fock->size;
	complex double A[s*s];
	uint ij = 0;
	for (uint i = 0; i < s; ++i)
		for (uint j = 0; j < s; ++j, ++ij)
			A[ij] = state_value_xy(fock->states[j], x[i], y[i]);
        double v = 1.0;
        for (uint i = 1; i <= s; ++i)
                v *= i;
	return permanent_z(A, s) / sqrt(v);
}

double fock_density(Fock *fock , double *x, double *y) {
	double v = fabs(fock_compute(fock, x, y));
	return v*v;
}

double fock_phase(Fock *fock , double *x, double *y) {
	return carg(fock_compute(fock, x, y));
}

/* === Way of handling sets of fock states === */
/* =========================================== */

void __fock_basis_alloc_recursive(FockBasis* set, Fock* running, uint n_particles, double e_cutoff) {
        if (n_particles == 0) {
                fock_basis_add(set, running);
                return;
        }
        for (uint i = 0; i <= running->states[n_particles]; ++i) {
                running->states[n_particles-1] = i;
                if (fock_energy(running) <= e_cutoff)
                        __fock_basis_alloc_recursive(set, running, n_particles-1, e_cutoff);
        }
        running->states[n_particles-1] = 0; /* Reset to the least energy state */
}

FockBasis fock_basis_alloc(uint n_particles, double e_cutoff) {
        FockBasis set = { 0 };
        Fock tmp = {
                .size = n_particles,
                .states = { 0 },
        };

        uint max = 0;
        while (state_energy(max) < e_cutoff)
                ++max;

        for (uint i = 0; i < max; ++i) {
                tmp.states[n_particles-1] = i;
                if (fock_energy(&tmp) <= e_cutoff)
                        __fock_basis_alloc_recursive(&set, &tmp, n_particles-1, e_cutoff);
        }
        return set;
}

FockBasis fock_basis_alloc_copy(FockBasis src) {
        FockBasis ret = {
                .states = util_malloc(sizeof(*src.states)*src.size),
                .size = src.size,
                .capacity = src.size,
        };
        memcpy(ret.states, src.states, sizeof(*src.states)*src.size);
        return ret;
}

void fock_basis_add(FockBasis* set, const Fock* fock) {
	++set->size;
        if (set->capacity < set->size) {
                set->capacity = MAX(set->size, 2*set->capacity);
	        set->states = util_realloc(set->states, sizeof(*set->states) * set->capacity);
        }
	memcpy(set->states + (set->size - 1), fock, sizeof(*fock));
}

void fock_basis_remove(FockBasis* basis, uint idx) {
        if (idx >= basis->size)
                return;
        --basis->size;
        memmove(&basis->states[idx], &basis->states[basis->size], sizeof(*basis->states));
}

void fock_basis_free(FockBasis* set) {
	free(set->states);
	set->size = 0;
        set->capacity = 0;
}

double fock_perturbation_energy_one(double g, Fock *state) {
        return g * fock_matrix_element_delta_potential(state, state);
}

double fock_perturbation_energy_two(double g, FockBasis left, Fock *state) {
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

double fock_perturbation_energy_three(double g, FockBasis set, Fock *state) {
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

uint fock_basis_find(FockBasis set, const Fock* state) {
	uint i = 0;
	for (; i < set.size; ++i)
		if (fock_equal(state, &set.states[i]))
			break;
	return i;
}

uint fock_operator_count(const Fock* fock, uint nr) {
	uint count = 0;
	for (uint i = 0; i < fock->size; ++i)
		if (fock->states[i] == nr)
			++count;
	return count;
}

double fock_operator_create(Fock* fock, uint nr) {
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

double fock_operator_annihilate(Fock* fock, uint nr) {
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

static void __fock_state_add_unnormalized(FockState fs[static 1], Fock fstate[static 1], complex double coeff) {
        const uint idx = fock_basis_find(fs->basis, fstate);
        if (idx != fs->basis.size) {
                fs->a[idx] += coeff;
                return;
        }

        /* There is fockstate in set, we have to add new one */
        uint old_cap = fs->basis.capacity;
        fock_basis_add(&fs->basis, fstate);
        if (old_cap != fs->basis.capacity)
                fs->a = util_realloc(fs->a, sizeof(*fs->a)*fs->basis.capacity);
        fs->a[idx] = coeff;
}

FockState fock_operator_psi(FockState fs, double x, double y) {
        /* Determining maximum annichilation operator */
        uint max = 0;
        for (uint i = 0; i < fs.basis.size; ++i)
                max = MAX(max, fs.basis.states[i].states[fs.basis.states[i].size-1]);
        ++max;

        FockState ret = {
                .a = util_malloc(sizeof(*ret.a)*16),
                .basis = {
                        .size = 0,
                        .capacity = 16,
                        .states = util_malloc(sizeof(*ret.basis.states)*16),
                },
        };


        Fock spstate;
        /* For every annichilation operator */
        for (uint i = 0; i < max; ++i) {
                const complex double vpsi = state_value(i, x, y);
                /* Apply the operator to every element of set */
                for (uint j = 0; j < fs.basis.size; ++j) {
                        complex double v  = fs.a[j];
                        if (FOCK_CONFIG_SET_WEIGHT_CUTOFF < cabs(v)) {
                                spstate = fs.basis.states[j];
                                v *= fock_operator_annihilate(&spstate, i) * vpsi;
                                if (v != 0.0)
                                        __fock_state_add_unnormalized(&ret, &spstate, v);
                        }
                }
        }

        fock_state_normalize(&ret);
        return ret;
}

void deigen_sol_free(DEigenSol *es) {
        free(es->U);
        free(es->val);
        es->size = 0;
}

void zeigen_sol_free(DEigenSol* es) {
        free(es->U);
        free(es->val);
        es->size = 0;
}

void fock_state_normalize(FockState* fs) {
        double tot = 0.0;
        for (uint i = 0; i < fs->basis.size; ++i) {
                // TODO: Check efficiency with creal() and cimag()
                double v = cabs(fs->a[i]);
                tot += v*v;
        }
        tot = sqrt(tot);
        for (uint i = 0; i < fs->basis.size; ++i)
                fs->a[i] /= tot;
}

void fock_state_compress(FockState* fs, double cutoff) {
        uint top = fs->basis.size;
        for (int i = fs->basis.size-1; i >= 0; --i)
                if (cabs(fs->a[i]) < cutoff) {
                        fock_basis_remove(&fs->basis, i);
                        fs->a[i] = fs->a[top];
                        --top;
                }

        /* TODO: Test the efficiency without reallocing */
        fs->basis.capacity = fs->basis.size;
        fs->basis.states = util_realloc(fs->basis.states, sizeof(*fs->basis.states)*fs->basis.size);
        fs->a = util_realloc(fs->basis.states, sizeof(*fs->a)*fs->basis.size);
        fock_state_normalize(fs);
}

FockState fock_state_alloc(FockBasis basis, double* a, double cutoff) {
        FockState ret = {
                .basis = fock_basis_alloc_copy(basis),
                .a = util_malloc(sizeof(*ret.a)*basis.size),
        };
        for (uint i = 0; i < basis.size; ++i)
                ret.a[i] = a[i];
        fock_state_compress(&ret, cutoff);
        return ret;
}

void fock_state_free(FockState* set) {
        fock_basis_free(&set->basis);
        free(set->a);
}

complex double fock_state_compute(FockState cset, double* x, double* y) {
	complex double ret = 0.0;
        //#pragma omp parallel for
	for (uint i = 0; i < cset.basis.size; ++i)
		if (FOCK_CONFIG_SET_WEIGHT_CUTOFF < cabs(cset.a[i]))
		        ret += cset.a[i] * fock_compute(&cset.basis.states[i], x, y);
	return ret;
}

double fock_matrix_element_kinetic(const Fock* left, const Fock* right) {
	if (!fock_equal(left, right))
		return 0.0;
	return fock_energy(left);
}

complex double fock_basis_compute(FockBasis set, double* x, double* y, const double* coeff) {
	complex double ret = 0.0;
        //#pragma omp parallel for
	for (uint i = 0; i < set.size; ++i)
		if (FOCK_CONFIG_SET_WEIGHT_CUTOFF < fabs(coeff[i]))
			ret += coeff[i] * fock_compute(&set.states[i], x, y);
	return ret;
}

double fock_basis_compute_density(FockBasis set, double* x, double* y, const double* coeff) {
	const double v = fabs(fock_basis_compute(set, x, y, coeff));
	return v*v;
}

double fock_stateset_compute_phase_weighted(FockBasis set, double* x, double* y, const double* coeff) {
	return carg(fock_basis_compute(set, x, y, coeff));
}

void fock_matrix_fill(FockBasis basis, double* T, double* V) {
        const uint s = basis.size;
        uint max = 0;
        for (uint i = 0; i < s; ++i)
                for (uint j = 0; j < basis.states[i].size; ++j)
                        max = MAX(basis.states[i].states[j], max);
        state_update_quick_functions(max);
        #pragma omp parallel for collapse(2)
        for (uint i = 0; i < s; ++i) {
                for (uint j = i; j < s; ++j) {
                        double t = fock_matrix_element_kinetic(&basis.states[j], &basis.states[i]);
                        double v = fock_matrix_element_delta_potential(&basis.states[j], &basis.states[i]);
                        V[i*s + j] = v;
                        V[j*s + i] = v;
                        T[i*s + j] = t;
                        T[j*s + i] = t;
                }
        }
}

DEigenSol fock_eigen_vt_solve(double g, uint n, double* T, double* V) {
        DEigenSol ret = {
                .size = n,
                .U = util_malloc(sizeof(*ret.U)*n*n),
                .val = util_malloc(sizeof(*ret.val)*n)
        };
        for (uint i = 0; i < n; ++i) {
                for (uint j = 0; j < n; ++j) {
                        ret.U[j*n+i] = T[j*n+i] + g*V[j*n+i];
                }
        }
	LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'L', ret.size, ret.U, ret.size, ret.val);
	return ret;
}

DEigenSol fock_eigen_solve(double g, FockBasis basis) {
        const uint s = basis.size;
        uint max = 0;
        for (uint i = 0; i < s; ++i)
                for (uint j = 0; j < basis.states[i].size; ++j)
                        max = MAX(basis.states[i].states[j], max);
        state_update_quick_functions(max);

	DEigenSol ret = {
		.size = s,
		.U = util_malloc(sizeof(*ret.U)*s*s),
		.val = util_malloc(sizeof(*ret.val)*s)
	};
        #pragma omp parallel for collapse(2)
	for (uint i = 0; i < s; ++i) {
		for (uint j = i; j < s; ++j) {
                        //printf("%3u %3u\n", i, j); fflush(stdout);
			const double t = fock_matrix_element_kinetic(&basis.states[j], &basis.states[i]);
			const double v = fock_matrix_element_delta_potential(&basis.states[j], &basis.states[i]);
			ret.U[i*s + j] = t + g*v;
		}
	}
	LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'L', s, ret.U, s, ret.val);
	return ret;
}

void fock_fprint(FILE* f, const Fock* fock) {
	for (uint i = 0; i < fock->size;) {
		state_fprint(f, fock->states[i]);
		if (++i != fock->size)
			putc(' ', f);
	}
}

void fock_basis_fprint(FILE* f, const FockBasis set) {
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

void fock_basis_fprint_weighted(FILE* f, const FockBasis set, const double* coeff, uint n) {
	size_t idxs[set.size];
	gsl_heapsort_index(idxs, coeff, set.size, sizeof(*coeff), __fock_square_compare);
	fprintf(f, " |ci|^2         state\n");
	for (uint i = 0; i < MIN(set.size, n); ++i) {
		fprintf(f, "%lf : ", coeff[idxs[i]]*coeff[idxs[i]]);
		fock_fprint(f, &set.states[idxs[i]]);
		putc('\n', f);
	}
}

void fock_coeffstateset_fprint(FILE* f, const FockState cset, uint n) {
	size_t idxs[cset.basis.size];
	gsl_heapsort_index(idxs, cset.a, cset.basis.size, sizeof(*cset.a), __fock_complex_square_compare);
	fprintf(f, " |ci|^2         state\n");
	for (uint i = 0; i < MIN(cset.basis.size, n); ++i) {
		fprintf(f, "%lf : ", cabs(cset.a[idxs[i]])*cabs(cset.a[idxs[i]]));
		fock_fprint(f, &cset.basis.states[idxs[i]]);
		putc('\n', f);
	}
}

complex double fock_state_aiaj(FockState fs, uint i, uint j) {
        FockBasis fb_copy = fock_basis_alloc_copy(fs.basis);
        double op_val[fb_copy.size];
        for (uint n = 0; n < fb_copy.size; ++n) {
                op_val[n]  = fock_operator_annihilate(fb_copy.states+n, j);
                op_val[n] *= fock_operator_create(fb_copy.states+n, i);
        }
        complex double ret = 0.0;
        for (uint n = 0; n < fs.basis.size; ++n) {
                if (op_val[n] == 0.0)
                        continue;

                for (uint m = 0; m < fs.basis.size; ++m) {
                        if (fock_equal(&fb_copy.states[n], &fs.basis.states[m])) {
                                ret += fs.a[n]*conj(fs.a[m])*op_val[n];
                                break;
                        }
                }
        }
        fock_basis_free(&fb_copy);
        return ret;
}

/* Procedures to recursivly genarate basis of fixed momentum fock states */

/* NOTE: One possible way of improvement is to check
 * if it is even possible to reached desired total momentum */

static void __basis_recursive(FockBasis* basis, int m_tot, int m_curr, double e_cut, double e_curr, Fock* fock, uint npart) {
        if (npart == 0) {
                if (m_tot == m_curr)
                        fock_basis_add(basis, fock);
                return;
        }

        for (uint i = 0; i <= fock->states[npart]; ++i) {
                fock->states[npart-1] = i;
                int m_next = m_curr + state_info(i).m;
                double e_next = e_curr + state_energy(i);
                if (e_curr <= e_cut)
                        __basis_recursive(basis, m_tot, m_next, e_cut, e_next, fock, npart-1);
        }
}

FockBasis fock_basis_alloc_fixed_angular_momentum(int m_tot, double e_cutoff, uint n_particles) {
        util_error(n_particles < 2, "Cannot generate states with < 2 particles.\n");

        uint max = 0;
        while (state_energy(max) < e_cutoff)
                ++max;

        FockBasis ret = { 0 };
        Fock fock = {
                .size = n_particles,
        };

        for (uint i = 0; i < max; ++i) {
                fock.states[n_particles-1] = i;
                int m_curr = state_info(i).m;
                double e_tot = state_energy(i);
                __basis_recursive(&ret, m_tot, m_curr, e_cutoff, e_tot, &fock, n_particles-1);
        }

        return ret;
}

double* fock_state_alloc_random(FockState fs, uint n) {
        uint npart = fs.basis.states[0].size;
        double* p = util_malloc(sizeof(*p)*n*2*npart);

        /* Draw first positions */
        // printf("# Allocating first density matrix.\n");
        Rho rho = rho_alloc(fs);
        // printf("# Drawing first positions\n");
        rho_rand(rho, FOCK_RHO_RAND_RELOAD, p, &p[npart]);
        for (uint i = 1; i < n; ++i)
                rho_rand(rho, FOCK_RHO_RAND_KEEP_ARR, &p[i*(2*npart)], &p[i*(2*npart)+npart]);
        rho_free(&rho);

        if (npart == 1)
                return p;

        Rho rhos[npart-1];
        FockState fss[npart];
        fss[0] = fs;
        // printf("# Drawing rest of positions.\n");
        for (uint i = 0; i < n; ++i) {
                fprintf(stderr, " [ %4u / %4u ]\n", i+1, n);
                for (uint j = 0; j < npart-1; ++j) {
                        //printf("%u", j+1);
                        fss[j+1] = fock_operator_psi(fss[j], p[i*(2*npart)+j], p[i*(2*npart)+npart+j]);
                        rhos[j] = rho_alloc(fss[j+1]);
                        rho_rand(rhos[j], FOCK_RHO_RAND_RELOAD, &p[i*(2*npart)+j+1], &p[i*(2*npart)+npart+j+1]);
                }
                //putchar('\n');


                for (uint j = 0; j <npart-1; ++j) {
                        fock_state_free(fss+j+1);
                        rho_free(rhos+j);
                }
        }

        return p;
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
