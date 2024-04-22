#include "fock.h"
#include "util.h"
#include <tgmath.h>

static uint _tot_num_planewave_states = 0;
static PlaneWaveState* _planewave_states;

static void __update_planewave_states(uint state_nr) {
	uint new_num = 0;
	int p = -1; // Maximal momentum
	for (uint i = 1;  state_nr >= new_num; ++p, i += 2)
		new_num = i*i;

	_planewave_states = util_realloc(_planewave_states, new_num * sizeof(PlaneWaveState));
	util_error(_planewave_states == NULL, "Failed to realloc\n");

	uint idxs[p+1];
	idxs[0] = 0;
	for (uint i = 1; i < p+1; ++i)
		idxs[i] = (i+i-1)*(i+i-1);
	for (int ix = -p; ix <= p; ++ix) for (int iy = -p; iy <= p; ++iy) {
		const int i = MAX(abs(ix), abs(iy));
		_planewave_states[idxs[i]].kx = ix;
		_planewave_states[idxs[i]].ky = iy;
		++idxs[i];
	}
	_tot_num_planewave_states = new_num;
}

double __planewave_energy(uint state_nr) {
	PlaneWaveState state = fock_planewave_state_index(state_nr);
	return 2.0*M_PI*M_PI*(state.kx*state.kx + state.ky*state.ky);
}

PlaneWaveState fock_planewave_state_index(uint state_nr) {
	if (_tot_num_planewave_states <= state_nr)
		__update_planewave_states(state_nr);
	return _planewave_states[state_nr];
}

uint fock_planewave_find_index(PlaneWaveState state) {
	uint max_k = MAX(abs(state.kx), abs(state.ky));
	uint i = (max_k == 0) ? 0 : (max_k+max_k-1)*(max_k+max_k-1);
	for (;;++i) {
		PlaneWaveState com = fock_planewave_state_index(i);
		if (com.kx == state.kx && com.ky == state.ky)
			break;
	}
	return i;
}

double complex __planewave_value(uint state_nr, double x, double y) {
	PlaneWaveState state = fock_planewave_state_index(state_nr);
	return exp(2*M_PI*I*(state.kx*x + state.ky*y));
}

Momentum fock_planewave_momentum(const FockState* fock) {
	util_error(fock->type != STATE_PLANEWAVE, "Cannot get momentum of state\n");
	Momentum r = { .kx = 0, .ky = 0 };
	for (uint i = 0; i < fock->size; ++i) {
		const Momentum m = fock_planewave_state_index(fock->states[i]);
		r.kx += m.kx;
		r.ky += m.ky;
	}
	return r;
}


double __fock_planewave_matrix_element_delta_potential(const FockState* left, const FockState* right) {
	const Momentum m_l = fock_planewave_momentum(left);
	const Momentum m_r = fock_planewave_momentum(right);
	if (m_l.kx != m_r.kx || m_l.ky != m_r.ky)
		return 0.0;

	FockState after1, after2, after3, after4;

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
					if (fock_equal(&after4, left))
						ret += v1*v2*v3*v4;

				}
			}
		}
	}

	return 0.5*ret;
}

void __fock_planewave_fprint(FILE* f, uint state_nr) {
	const PlaneWaveState pw = fock_planewave_state_index(state_nr);
	fprintf(f, "( %2d, %2d )", pw.kx, pw.ky);
}

Vec2z __fock_planewave_grad(uint state_nr, double x, double y, complex double* gx, complex double* gy) {
        PlaneWaveState ws = fock_planewave_state_index(state_nr);
        const complex double val = __planewave_value(state_nr, x, y);
        Vec2z ret = {
                x = I*ws.kx*val,
                y = I*ws.ky*val
        };
        return ret;
}

