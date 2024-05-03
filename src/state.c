#include <string.h>
#include <tgmath.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_heapsort.h>

#include "fock.h"
#include "state.h"
#include "util.h"



/* Static elements used for quick ordering of single particle states */

static uint _bessel_num_zeros = 0;
static uint _bessel_tot_num_states = 0;
static StateInfo* _bessel_states;

/* Arrays used for quick integrations of bessel fucntions,
 * those are updated after calling functions that give values */

static uint _bessel_num_qf = 0;
static QuickFunction *_bessel_qf_array;

StateInfo state_info_direct(int m, uint zero_num) {
        util_error(zero_num == 0, "Cannot have 0-th 'zero of bessel funtion'.\n");
        StateInfo ret = {
                .m = m,
                .zero_num = zero_num,
                .zero_val = gsl_sf_bessel_zero_Jnu(abs(m), zero_num),
                .normalization_constant = 1/(fabs(gsl_sf_bessel_Jn(abs(m)+1, ret.zero_val))*sqrt(M_PI))
        };
        return ret;
}

int fock_angular_momentum(const Fock* state) {
        int m_tot = 0;
        for (uint i = 0; i < state->size; ++i)
                m_tot += state_info(state->states[i]).m;
        return m_tot;
}

static int __bessel_state_compare(const void* l, const void* r) {
        StateInfo* ls = (StateInfo*) l;
        StateInfo* rs = (StateInfo*) r;
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
                        _bessel_states[ixy] = state_info_direct(mom, n_zero);
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

StateInfo state_info(uint state_nr) {
        __bessel_ensure_enough(state_nr);
	return _bessel_states[state_nr];
}

double state_energy(uint n) {
        __bessel_ensure_enough(n);
        const double k = _bessel_states[n].zero_val;
	return 0.5*k*k;
}

void state_update_quick_functions(uint max) {
        if (max < _bessel_num_qf)
                return;
        __bessel_ensure_enough(max+9);
        _bessel_qf_array = util_realloc(_bessel_qf_array, sizeof(*_bessel_qf_array)*(max+10));
        for (uint i = _bessel_num_qf; i < max+10; ++i) {
                StateInfo bws = _bessel_states[i];
                for (uint j = 0; j < QF_ARR_LEN; ++j) {
                        double r = j * QF_DR;
                        double val = gsl_sf_bessel_Jn(bws.m, bws.zero_val*r);
                        _bessel_qf_array[i].arr[j] = bws.normalization_constant * val;
                }
        }
        _bessel_num_qf = max + 10;
}

double state_value_r(uint n, double r) {
        state_update_quick_functions(n);
        return qf_get_value(_bessel_qf_array+n, r);
}

double complex state_value(uint n, double r, double phi) {
        double rval = state_value_r(n, r);
        StateInfo si = _bessel_states[n];
        return rval * cexp(I*si.m*phi);
}

double complex state_value_xy(uint state_nr, double x, double y) {
        const double r = sqrt(x*x + y*y);
        const double phi = atan2(y, x);
        return state_value(state_nr, r, phi);
}

const QuickFunction* state_get_qf(uint nr) {
        state_update_quick_functions(nr);
        return &_bessel_qf_array[nr];
}

static double __bessel_quadruple_integral(uint s1, uint s2, uint s3, uint s4) {
        QuickFunction iqf = *state_get_qf(s1);
        const QuickFunction *qf1 = state_get_qf(s2);
        const QuickFunction *qf2 = state_get_qf(s3);
        const QuickFunction *qf3 = state_get_qf(s4);
        for (uint i = 0; i < QF_ARR_LEN; ++i)
                iqf.arr[i] *= qf1->arr[i] * qf2->arr[i] * qf3->arr[i];
        return qf_integrate(&iqf);
}

double fock_matrix_element_delta_potential(const Fock* left, const Fock* right) {
        int ml = fock_angular_momentum(left);
        int mr = fock_angular_momentum(right);
        if (ml != mr)
                return 0.0;

	Fock after1;
	Fock after2;
	Fock after3;
	Fock after4;

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

        /* Coefficient comes from 1/2 before sum in interaction term
         * and 2*M_PI from integral over angle */
	return M_PI*ret;
}

void state_fprint(FILE* f, uint state_nr) {
        StateInfo bws = state_info(state_nr);
        fprintf(f, "(%d, m:%3d)", bws.zero_num, bws.m);
}

uint state_find_index(StateInfo state) {
        for (uint i = 0; i < 10000;++i) {
                StateInfo what = state_info(i);
                if (what.m == state.m && what.zero_num == state.zero_num)
                        return i;
        }
        util_error(true, "Cannot find singleparticle state index.\n");
        return 0;
}
