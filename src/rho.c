#include <tgmath.h>
#include <string.h>

#include <cblas.h>
#include <lapacke.h>

#include "rho.h"
#include "util.h"

Rho rho_alloc(FockState fs) {
        Rho ret = { 0 };
        for (uint i = 0; i < fs.basis.size; ++i)
                for (uint j = 0; j < fs.basis.states[i].size; ++j)
                        if (!state_basis_contains(&ret.basis, fs.basis.states[i].states[j]))
                                state_basis_add(&ret.basis, fs.basis.states[i].states[j]);


        const uint s = ret.basis.size;
        ret.M = util_malloc(sizeof(*ret.M)*s*s);

        // TODO: Faster way of computing aiaj elements of a rho matrix
        for (uint i = 0; i < s; ++i) {
                const uint si = ret.basis.states[i];
                for (uint j = i; j < s; ++j) {
                        const uint sj = ret.basis.states[j];
                        complex double v = fock_state_aiaj(fs, sj, si);
                        ret.M[i*s + j] = conj(v);
                        ret.M[j*s + i] = v;
                }
        }

        return ret;
}

void rho_free(Rho* rho) {
        state_basis_free(&rho->basis);
        free(rho->M);
}

complex double rho_compute(Rho rho, double r, double phi, double rp, double phip) {
        // TODO: Can we use the fact that rho is hermitian?
        
        const uint s = rho.basis.size;
        complex double vec[s];
        for (uint i = 0; i < s; ++i)
                vec[i] = state_value(rho.basis.states[i], r, phi);
        cblas_ztrmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, s, rho.M, s, vec, 1);
        complex double ret = 0.0;
        for (uint i = 0; i < s; ++i)
                ret += vec[i] * conj(state_value(rho.basis.states[i], rp, phip));

        return ret;
}

complex double rho_compute_xy(Rho rho, double x, double y, double xp, double yp) {
        double r = sqrt(x*x + y*y);
        double phi = atan2(y, x);
        double rp = sqrt(xp*xp + yp*yp);
        double phip = atan2(yp, xp);

        return rho_compute(rho, r, phi, rp, phip);
}

double rho_compute_diag(Rho rho, double r, double phi) {
        const uint s = rho.basis.size;
        complex double vec[s];
        for (uint i = 0; i < s; ++i)
                vec[i] = state_value(rho.basis.states[i], r, phi);
        double ret = 0.0;
        for (uint i = 0; i < s; ++i) {
                complex double sum = 0.0;
                for (uint j = i+1; j < s; ++j) {
                        sum += rho.M[j*s + i]*vec[j];
                }
                ret += 2*creal(conj(vec[i])*sum);
                double v = cabs(vec[i]);
                ret += v*v*creal(rho.M[i*s + i]);
        }

        return ret;
}

double rho_compute_diag_xy(Rho rho, double x, double y) {
        double r = sqrt(x*x + y*y);
        double phi = atan2(y, x);
        return rho_compute_diag(rho, r, phi);
}

DEigenSol rho_diag(Rho* rho) {
        const uint s = rho->basis.size;
        DEigenSol ret = {
                .U = util_malloc(sizeof(*ret.U)*s*s),
                .val = util_malloc(sizeof(*ret.val)*s),
                .size = s,
        };
        memcpy(ret.U, rho->M, sizeof(*ret.U)*s*s);
        LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'L', s, ret.U, s, ret.val);
        return ret;
}

static uint __bsearch_double(double key, double* arr, uint n) {
        for (uint i = 0; i < n; ++i)
                if (key <= arr[i])
                        return i;

        return n-1;
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

void rho_rand(Rho rho, RhoRandOpt opt, double *x, double *y) {
        const double dr   = 1.0 / _G1_NR;
        const double dphi = 2.0 * M_PI / _G1_NPHI;

        /* Only if proper option is set, we recalculate the integral */
        if (opt == FOCK_RHO_RAND_RELOAD) {
                double rsum = 0.0;
                for (uint i = 0; i < _G1_NR; ++i) {
                        const double r = (i+0.5)*dr; 
                        const double surf = 0.5*dphi*dr*dr*(2*i+1);
                        double sum = 0.0;
                        for (uint j = 0; j < _G1_NPHI; ++j) {
                                const double phi = (j+0.5)*dphi; 
                                sum += surf*rho_compute_diag(rho, r, phi);
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

        *x = r*cos(phi);
        *y = r*sin(phi);
}
