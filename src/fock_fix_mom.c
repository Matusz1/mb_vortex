#include <stdio.h>
#include <string.h>
#include "util.h"
#include "fock.h"

#define ABS(x) ((x) < 0 ? (-(x)) : (x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))

typedef struct {
	uint* indices;    /* Used to index vectores */
	uint n_particles;
	uint num;         /* There is num*n_particles indices */
	uint capacity;    /* num < capacity */
} VecCollection;

bool are_equal(uint* v1, uint* v2, uint n) {
	return memcmp(v1, v2, n*sizeof(*v1)) == 0;
}

void collection_add(VecCollection* coll, uint* ind) {
	if (coll->num == coll->capacity) {
		coll->capacity = MAX(16, coll->capacity*2);
		coll->indices = util_realloc(coll->indices, sizeof(*coll->indices) * coll->capacity * coll->n_particles);
	}
	memcpy(coll->indices+coll->num*coll->n_particles, ind, sizeof(*ind)*coll->n_particles);
	++coll->num;
}

void collection_add_unique(VecCollection* coll, uint* ind) {
	for (uint i = 0; i < coll->num; ++i)
		if (are_equal(ind, coll->indices + i*coll->n_particles, coll->n_particles))
			return;
	collection_add(coll, ind);
}

bool is_reachable(Momentum dest, Momentum cp, int max) {
	return (ABS(dest.kx-cp.kx) <= max) && (ABS(dest.ky-cp.ky) <= max);
}

void generate_collection_recursive(VecCollection* coll, Momentum dest, Momentum cp, int e_cut, int e_tot, int max_d, uint* ind, uint n_part) {
	if (e_tot > e_cut || !is_reachable(dest, cp, max_d*n_part))
		return;
	if (n_part == 0) {
		if (dest.kx == cp.kx && dest.ky == cp.ky)
			collection_add_unique(coll, ind);
	} else {
		for (uint i = 0; i <= ind[n_part]; ++i) {
			ind[n_part-1] = i;
			Momentum mv = fock_planewave_state_index(i);
			uint next_ener = mv.kx*mv.kx + mv.ky*mv.ky + e_tot;
			mv.kx += cp.kx;
			mv.ky += cp.ky;
			generate_collection_recursive(coll, dest, mv, e_cut, next_ener, max_d, ind, n_part-1);
		}
	}
}

VecCollection generate_collection(Momentum dst, int e_cutoff, uint max_state, uint n_particles) {
        uint size = (max_state+max_state+1)*(max_state+max_state+1);
	VecCollection collection = {
		.num = 0,
		.capacity = 16,
		.n_particles = n_particles,
		.indices = util_malloc(sizeof(*collection.indices) * 16*n_particles)
	};
	uint curr[n_particles];
	memset(curr, 0, n_particles*sizeof(*curr));
	for (uint i = 0; i < size; ++i) {
		curr[n_particles-1] = i;
		Momentum cp = fock_planewave_state_index(i);
		int e_total = cp.kx*cp.kx + cp.ky*cp.ky;
		generate_collection_recursive(&collection, dst, cp, e_cutoff, e_total, max_state, curr, n_particles-1);
	}

	return collection;
}

FockStateSet fock_stateset_alloc_same_momentum(int kx, int ky, int e_cutoff, uint max_state, uint n_particles) {
        Momentum dst = { .kx = kx, .ky = ky };
        VecCollection col = generate_collection(dst, e_cutoff, max_state, n_particles);
        FockStateSet set = {
                .size = col.num,
                .states = util_malloc(sizeof(*set.states) * col.num)
        };
        for (uint i = 0; i < col.num; ++i) {
                memcpy(set.states[i].states, col.indices + i*col.n_particles, col.n_particles*sizeof(*col.indices));
                set.states[i].size = n_particles;
                set.states[i].type = STATE_PLANEWAVE;
        }
        util_free(col.indices);
        return set;
}
