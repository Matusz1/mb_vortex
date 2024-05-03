#include <tgmath.h>
#include <string.h>

#include <lapacke.h>
#include <cblas.h>

#include "gpe_radial.h"
#include "fock.h"
#include "state.h"
#include "util.h"

double __gpe_solve_calc_integral(const QuickFunction* chi2, uint nr_1, uint nr_2) {
        const QuickFunction qf1 = *state_get_qf(nr_1);
        const QuickFunction qf2 = *state_get_qf(nr_2);
        //util_error(qf1 == NULL, "NULL pointer from state_get_qf(nr_1)\n");
        //util_error(qf2 == NULL, "NULL pointer from state_get_qf(nr_2)\n");
        QuickFunction ifun;
        for (uint i = 0; i < QF_ARR_LEN; ++i) {
                const double v1 = qf1.arr[i];
                const double v2 = qf2.arr[i];
                ifun.arr[i] = chi2->arr[i] * v1 * v2;
        }
        return qf_integrate(&ifun);
}

BesselLinear gpe_solve_m1(double gN, double* mu, double* Epp) {
        const uint s = BESSEL_LINEAR_COMB_MAX_SIZE;
        BesselLinear ret = {
                .size = s,
        };
        for (uint i = 0; i < s; ++i) {
                ret.w[i] = 0.0;
                StateInfo si = state_info_direct(1, i+1);
                ret.nrs[i] = state_find_index(si);
        }
        ret.w[0] = 1.0;

        QuickFunction qf_chi2 = *state_get_qf(ret.nrs[0]);
        for (uint i = 0; i < QF_ARR_LEN; ++i)
                qf_chi2.arr[i] *= qf_chi2.arr[i]; 


        printf("# GPE selfconsistent diagonalization starting ...\n");

        double H[s*s];
        double diff = 1.0;
        *mu = 1.0;
        double mis[s];

        uint counter = 1;
        while (diff > 0.0000001 && counter <= 10000) {
                /* Fill Hamiltonian */
                for (uint n = 0; n < s; ++n) {
                        for (uint m = n; m < s; ++m) {
                                double itgr = __gpe_solve_calc_integral(&qf_chi2, ret.nrs[n], ret.nrs[m]);
                                H[n*s+m] = 2*M_PI*gN*itgr;
                        }
                        H[n*s+n] += state_energy(ret.nrs[n]);
                }

                LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', s, H, s, mis);

                diff = fabs((mis[0]-*mu)/(*mu));
                *mu = mis[0];
                for (uint i = 0; i < s; ++i)
                        ret.w[i] = H[i*s];

                /* Determine Chi */
                for (uint i = 0; i < QF_ARR_LEN; ++i) {
                        double val = 0.0;
                        for (uint j = 0; j < s; ++j) {
                                const QuickFunction *qf = state_get_qf(ret.nrs[j]);
                                double ai = ret.w[j];
                                val += ai * qf->arr[i];
                        }
                        qf_chi2.arr[i] = 0.5*val*val + 0.5*qf_chi2.arr[i];
                }

                printf("# %5u:\tmu = %lf\tdiff = %e\n", counter, *mu, diff);
                ++counter;
        }

        printf("# GPE CONVERGED, FINISIHING...\n");
        double frac[10] = {
                0.5, 0.6, 0.7, 0.8, 0.9,
                0.9, 0.9, 1.0, 1.0, 1.0,
        };

        for (uint i = 0; i < 10; ++i) {
                /* Fill Hamiltonian */
                for (uint n = 0; n < s; ++n) {
                        for (uint m = n; m < s; ++m) {
                                double itgr = __gpe_solve_calc_integral(&qf_chi2, ret.nrs[n], ret.nrs[m]);
                                H[n*s+m] = 2*M_PI*gN*itgr;
                        }
                        H[n*s+n] += state_energy(ret.nrs[n]);
                }

                LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', s, H, s, mis);

                *mu = mis[0];
                for (uint j = 0; j < s; ++j)
                        ret.w[j] = H[j*s];

                /* Determine Chi */
                for (uint j = 0; j < QF_ARR_LEN; ++j) {
                        double val = 0.0;
                        for (uint k = 0; k < s; ++k) {
                                const QuickFunction *qf = state_get_qf(ret.nrs[k]);
                                double ai = ret.w[k];
                                val += ai * qf->arr[j];
                        }
                        qf_chi2.arr[i] = frac[i]*val*val + (1.0-frac[i])*qf_chi2.arr[i];
                }

                printf("# %5u:\tmu = %lf\n", i, *mu);
        }

        for (uint i = 0; i < QF_ARR_LEN; ++i)
                qf_chi2.arr[i] *= qf_chi2.arr[i];

        double eshift = M_PI*gN*qf_integrate(&qf_chi2); 
        *Epp = *mu - eshift;
        printf("# GPE: SELF-CONSISTANT DONE\n");
        printf("# GPE: mu = %lf,    Epp = %lf\n", *mu, *Epp);
        
        if (ret.w[0] > 0)
                return ret;

        for (uint i = 0; i < s; ++i)
                ret.w[i] *= -1;

        return ret;
}

QuickFunction gpe_bes_lin_get_qf(const BesselLinear* bl) {
        QuickFunction ret = { 0 };
        for (uint i = 0; i < bl->size; ++i) {
                const QuickFunction *qf = state_get_qf(bl->nrs[i]);
                const double ai = bl->w[i];
                for (uint j = 0; j < QF_ARR_LEN; ++j)
                        ret.arr[j] += ai * qf->arr[j];
        }
        return ret;
}

static double __bdg_hamiltonian_elem(const QuickFunction* chi2, UVBesselState l, UVBesselState r) {
        const double integral = 2 * M_PI * __gpe_solve_calc_integral(chi2, l.nr, r.nr);
        if (l.isUp && r.isUp)
                return 2 * integral;
        if ((!l.isUp) && (!r.isUp))
                return -2 * integral;
        if ((!l.isUp) && r.isUp)
                return -integral;
        return integral;
}


UVSol bdg_solve(uint half_size, int m, double gN, double omega) {
        const uint size = 2 * half_size;
        double mu, Epp;
        BesselLinear blin = gpe_solve_m1(gN, &mu, &Epp);
        QuickFunction chi2 = gpe_bes_lin_get_qf(&blin);
        for (uint i = 0; i < QF_ARR_LEN; ++i) {
                chi2.arr[i] *= chi2.arr[i];
        }

        /* Define the states, fill non-iteraction part of Hamiltonian */
        double H[size*size];
        memset(H, 0, sizeof(*H)*size*size);
        UVSol sol = {
                .size = size,
                .types = util_malloc(sizeof(*sol.types)*size),
                .evals = util_malloc(sizeof(*sol.evals)*size),
                .evecs = util_malloc(sizeof(*sol.evecs)*size*size),
                .plus_family = util_malloc(sizeof(*sol.plus_family)*size),
        };
        for (uint i = 0; i < half_size; ++i) {
                StateInfo usi = state_info_direct(m+2, i+1);
                StateInfo lsi = state_info_direct(m, i+1);
                uint uidx = state_find_index(usi);
                uint lidx = state_find_index(lsi);

                sol.types[i].isUp = true;
                sol.types[i].nr = uidx;
                sol.types[i+half_size].isUp = false;
                sol.types[i+half_size].nr = lidx;

                H[i*size + i] = state_energy(uidx) - omega*(m+2) - mu;
                H[(i+half_size)*size + (i+half_size)] = -state_energy(lidx) + omega*m + mu;
        }

        /* Fill Interaction part */
        for (uint i = 0; i < size; ++i) {
                for (uint j = 0; j < size; ++j) {
                        H[i*size + j] += gN * __bdg_hamiltonian_elem(&chi2, sol.types[i], sol.types[j]);
                }
        }

        /* Find solution */
        double wr[size];
        double wi[size];
        double vr[size*size];
        LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', size, H, size, wr, wi, NULL, size, vr, size);

        /* Save the solution */
        for (uint i = 0; i < size; ++i) {
                if (fabs(wi[i]) == 0.0) {
                        sol.evals[i] = wr[i];
                        for (uint j = 0; j < size; ++j)
                                sol.evecs[j*size + i] = vr[j*size + i];
                } else {
                        sol.evals[i] = wr[i] + I*wi[i];
                        sol.evals[i+1] = wr[i+1] - I*wi[i+1];
                        for (uint j = 0; j < size; ++j) {
                                sol.evecs[j*size + i] = vr[j*size + i] + I*vr[j*size + i+1];
                                sol.evecs[j*size + i+1] = vr[j*size + i] - I*vr[j*size + i+1];
                        }
                        ++i;
                }
        }

        /* Determine family type */
        for (uint i = 0; i < size; ++i) {
                double pf = 0.0;
                double mf = 0.0;
                for (uint j = 0; j < size; ++j) {
                        double v = cabs(sol.evecs[j*size + i]);
                        v *= v;
                        if (j < half_size)
                                pf += v;
                        else
                                mf += v;
                }
                double pmi = pf - mf;
                if (pmi > 0)
                        sol.plus_family[i] = true;
                else
                        sol.plus_family[i] = false;
        }

        return sol;
}

void uvsol_free(UVSol* uvsol) {
        free(uvsol->types);
        free(uvsol->evals);
        free(uvsol->evecs);
}
