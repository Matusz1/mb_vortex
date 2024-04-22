#include "fock.h"

#include <string.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_heapsort.h>

#include "util.h"

#include "fock.h"

/* Static elements used for quick ordering of single particle states */

static uint _bessel_num_zeros = 0;
static uint _bessel_tot_num_states = 0;
static BesselWaveState* _bessel_states;

/* Arrays used for quick integrations of bessel fucntions,
 * those are updated after calling functions that give values */

static uint _bessel_num_qf = 0;
static QuickFunction *_bessel_qf_array;

BesselWaveState fock_bessel_state_direct(int m, uint zero_num) {
        util_error(zero_num == 0, "Cannot have 0-th 'zero of bessel funtion'.\n");
        BesselWaveState ret = {
                .m = m,
                .zero_num = zero_num,
                .zero_val = gsl_sf_bessel_zero_Jnu(abs(m), zero_num),
                .normalization_constant = 1/(fabs(gsl_sf_bessel_Jn(abs(m)+1, ret.zero_val))*sqrt(M_PI))
        };
        return ret;
}

int fock_bessel_angular_momentum(const FockState* state) {
        int m_tot = 0;
        for (uint i = 0; i < state->size; ++i)
                m_tot += fock_bessel_state_index(state->states[i]).m;
        return m_tot;
}

static int __bessel_state_compare(const void* l, const void* r) {
        BesselWaveState* ls = (BesselWaveState*) l;
        BesselWaveState* rs = (BesselWaveState*) r;
        if (ls->zero_val == rs->zero_val) {
                if (ls->m < rs->m)
                        return -1;
                if (ls->m == rs->m)
                        return 0;
        }
        if (ls->zero_val < rs->zero_val)
                return -1;
        return 1;
}

static void __bessel_realloc_states(uint nnz) {
        if (nnz == 0)
                return;
        int nx = nnz;
        int ny = 3*nnz;
	const uint e_size = nx * ny;
	double tmp_zero_pos[e_size];
	uint ixy = 0;
	for (uint i = 0; i < ny; ++i) {
		for (uint j = 0; j < nx; ++j) {
			tmp_zero_pos[ixy] = gsl_sf_bessel_zero_Jnu(i, j+1);
			++ixy;
		}
	}
	/* Bottom left has to smaller than top right */
        const double max_zero = tmp_zero_pos[nx-1];
	if (nx > 1)
		while (tmp_zero_pos[(ny-1)*nx] > max_zero)
			--ny;

        const uint matrix_size = 2*(ny-1)*nx + nx;
	_bessel_tot_num_states = matrix_size;
	_bessel_states = util_realloc(_bessel_states, sizeof(*_bessel_states) * _bessel_tot_num_states);
	_bessel_num_zeros = nx;

	ixy = 0;
	for (int mom = -ny+1; mom <= ny-1; ++mom) {
		for (uint n_zero = 1; n_zero <= nx; ++n_zero) {
                        _bessel_states[ixy] = fock_bessel_state_direct(mom, n_zero);
                        if (_bessel_states[ixy].zero_val > max_zero)
                                --_bessel_tot_num_states;
			++ixy;
		}
	}
        /* We have to sort the states using energy and momentum */
        gsl_heapsort(_bessel_states, matrix_size, sizeof(*_bessel_states), __bessel_state_compare);
}

static inline void __bessel_ensure_enough(uint state_nr) {
	while (_bessel_tot_num_states <= state_nr)
		__bessel_realloc_states(_bessel_num_zeros+1);
}

BesselWaveState fock_bessel_state_index(uint state_nr) {
        __bessel_ensure_enough(state_nr);
	return _bessel_states[state_nr];
}

double __bessel_energy(uint state_nr) {
        __bessel_ensure_enough(state_nr);
        const double k = _bessel_states[state_nr].zero_val;
	return 0.5*k*k;
}

static void __bessel_update_arrays(uint max) {
        if (max < _bessel_num_qf)
                return;
        __bessel_ensure_enough(max+9);
        _bessel_qf_array = util_realloc(_bessel_qf_array, sizeof(*_bessel_qf_array)*(max+10));
        for (uint i = _bessel_num_qf; i < max+10; ++i) {
                BesselWaveState bws = _bessel_states[i];
                for (uint j = 0; j < QF_ARR_LEN; ++j) {
                        double r = j * QF_DR;
                        double val = gsl_sf_bessel_Jn(bws.m, bws.zero_val*r);
                        _bessel_qf_array[i].arr[j] = bws.normalization_constant * val;
                }
        }
        _bessel_num_qf = max + 10;
}

double bessel_value_radial(uint state_nr, double r) {
        __bessel_update_arrays(state_nr);
        return qf_get_value(_bessel_qf_array+state_nr, r);
}

double complex bessel_value(uint state_nr, double r, double phi) {
        double rval = bessel_value_radial(state_nr, r);
        BesselWaveState bws = _bessel_states[state_nr];
        return rval * cexp(I*bws.m*phi);
}

double complex __bessel_value_xy(uint state_nr, double x, double y) {
        const double r = sqrt(x*x + y*y);
        const double phi = atan2(y, x);
        return bessel_value(state_nr, r, phi);
}

const QuickFunction* bessel_get_qf(uint nr) {
        __bessel_update_arrays(nr);
        return &_bessel_qf_array[nr];
}

/* params - uint-s with indexes of bessel functions */
static double __bessel_value_bessel_quadruple(double r, void* params) {
        double ret = r;
        for (uint i = 0; i < 4; ++i)
                ret *= bessel_value_radial(((uint*)params)[i], r);
        return ret;
}

static double __bessel_quadruple_integral(uint s1, uint s2, uint s3, uint s4) {
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
        gsl_function func;
        uint params[4] = { s1, s2, s3, s4 };
        func.function = __bessel_value_bessel_quadruple;
        func.params = params;
        double result, error;
        gsl_integration_qags(&func, 0, 1, 0, 1e-7, 1000, w, &result, &error);
        gsl_integration_workspace_free(w);
        return result;
}

double __bessel_matrix_element_delta_potential(const FockState* left, const FockState* right) {
        int ml = fock_bessel_angular_momentum(left);
        int mr = fock_bessel_angular_momentum(right);
        if (ml != mr)
                return 0.0;

	FockState after1;
	FockState after2;
	FockState after3;
	FockState after4;

	double ret = 0.0;
	/* This is a bit tricky because we have to iterate over every UNIQUE state */
	uint prev1 = right->states[0];
	uint prev2 = right->states[0];
	uint prev3 = left->states[0];
	uint prev4 = left->states[0];
	for (uint i = 0; i < right->size; ++i) {
		if (i != 0 && prev1 == right->states[i])
			continue;
                after1 = *right;
		double v1 = fock_operator_annihilate(&after1, right->states[i]);
		prev1 = right->states[i];
		for (uint j = 0; j < after1.size; ++j) {
			if (j != 0 && prev2 == after1.states[j])
				continue;
                        after2 = after1;
			double v2 = fock_operator_annihilate(&after2, after1.states[j]);
			prev2 = after1.states[j];
			for (uint q = 0; q < left->size; ++q) {
				if (q != 0 && prev3 == left->states[q])
					continue;
                                after3 = after2;
				double v3 = fock_operator_create(&after3, left->states[q]);
				prev3 = left->states[q];
				for (uint p = 0; p < left->size; ++p) {
					if (p != 0 && prev4 == left->states[p])
						continue;
                                        after4 = after3;
					double v4 = fock_operator_create(&after4, left->states[p]);
					prev4 = left->states[p];
					if (fock_equal(&after4, left)) {
                                                double ibes = __bessel_quadruple_integral(prev1, prev2, prev3, prev4);
						ret += v1*v2*v3*v4*ibes;
                                        }

				}
			}
		}
	}

	return M_PI*ret;
}

void __fock_bessel_fprint(FILE* f, uint state_nr) {
        BesselWaveState bws = fock_bessel_state_index(state_nr);
        fprintf(f, "(%d, m:%3d)", bws.zero_num, bws.m);
}

Vec2z __fock_bessel_grad(uint state_nr, double x, double y) {
        BesselWaveState ws = fock_bessel_state_index(state_nr);
        const double r = sqrt(x*x+y*y);
        if (r > 1.0)
                return (Vec2z) { 0 };
        const double phi = 2*atan(y/(x+r));
        complex double dr = ws.normalization_constant * ws.zero_val * 0.5;
        dr *= (gsl_sf_bessel_Jn(ws.m-1, ws.zero_val*r) - gsl_sf_bessel_Jn(ws.m+1, ws.zero_val*r)) * cexp(I*ws.m*phi);
        complex double dphi = I*ws.m*ws.normalization_constant*gsl_sf_bessel_Jn(ws.m, ws.zero_val*r) * cexp(I*ws.m*phi) / r;
        Vec2z ret = {
                .x = cos(phi)*dr - sin(phi)*dphi,
                .y = sin(phi)*dr + cos(phi)*dphi
        };
        return ret;
}

uint fock_bessel_find_index(BesselWaveState state) {
        for (uint i = 0; i < 10000;++i) {
                BesselWaveState what = fock_bessel_state_index(i);
                if (what.m == state.m && what.zero_num == state.zero_num)
                        return i;
        }
        util_error(true, "Cannot find singleparticle state index.\n");
        return 0;
}

/* === Routines for generating set of states with fixed momentum === */

typedef struct {
	uint* indices;    /* Used to index vectores */
	uint n_particles;
	uint num;         /* There is num*n_particles indices */
	uint capacity;    /* num < capacity */
} VecCollection;

static inline bool are_equal(uint* v1, uint* v2, uint n) {
	return memcmp(v1, v2, n*sizeof(*v1)) == 0;
}

static void collection_add(VecCollection* coll, uint* ind) {
	if (coll->num == coll->capacity) {
		coll->capacity = MAX(16, coll->capacity*2);
		coll->indices = util_realloc(coll->indices, sizeof(*coll->indices) * coll->capacity * coll->n_particles);
	}
	memcpy(coll->indices+coll->num*coll->n_particles, ind, sizeof(*ind)*coll->n_particles);
	++coll->num;
}

static inline bool __gen_is_reachable(int dest, int cp, int max) {
	return (abs(dest-cp) <= max);
}

static void generate_collection_recursive(VecCollection* coll, int m_dest, int m_curr, double e_cut, double e_tot, int max_d, uint* ind, uint n_part) {
	if (!__gen_is_reachable(m_dest, m_curr, max_d*n_part))
		return;
	if (n_part == 0) {
		if (m_dest == m_curr)
			collection_add(coll, ind);
	} else {
		for (uint i = 0; i <= ind[n_part]; ++i) {
			ind[n_part-1] = i;
                        BesselWaveState bws = fock_bessel_state_index(i);
                        double next_ener = __bessel_energy(i) + e_tot;
                        if (next_ener <= e_cut)
			        generate_collection_recursive(coll, m_dest, m_curr+bws.m, e_cut, next_ener, max_d, ind, n_part-1);
		}
	}
}

static VecCollection generate_collection(int m_dest, double e_cutoff, uint max_state, uint n_particles) {
        /* Find maximum single particle momentum */
        int m_max = 0;
        for (uint i = 0; i < max_state; ++i) {
                int m_curr = abs(fock_bessel_state_index(i).m);
                if (m_curr > m_max)
                        m_max = m_curr;
        }
	VecCollection collection = {
		.num = 0,
		.capacity = 16,
		.n_particles = n_particles,
		.indices = util_malloc(sizeof(*collection.indices) * 16*n_particles)
	};
	uint curr[n_particles];
	for (uint i = 0; i < max_state; ++i) {
		curr[n_particles-1] = i;
		int m_curr = fock_bessel_state_index(i).m;
                double e_total = __bessel_energy(i);
		generate_collection_recursive(&collection, m_dest, m_curr, e_cutoff, e_total, m_max, curr, n_particles-1);
	}

	return collection;
}

FockStateSet fock_stateset_alloc_same_angular_momentum(int m_tot, double e_cutoff, uint max_index, uint n_particles) {
        util_error(n_particles < 2, "Cannot generate states with < 2 particles.\n");
        VecCollection col = generate_collection(m_tot, e_cutoff, max_index, n_particles);
        FockStateSet set = {
                .size = col.num,
                .states = util_malloc(sizeof(*set.states) * col.num)
        };
        for (uint i = 0; i < col.num; ++i) {
                memcpy(set.states[i].states, col.indices + i*col.n_particles, col.n_particles*sizeof(*col.indices));
                set.states[i].size = n_particles;
                set.states[i].type = STATE_BESSEL;
        }
        util_free(col.indices);
        return set;

}
