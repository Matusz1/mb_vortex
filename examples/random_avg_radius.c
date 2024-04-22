#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include "src/fock.h"
#include <gsl/gsl_sf_bessel.h>

#include <rgc_memory.h>
#include <rgc_io.h>

int main(int argc, char **argv) {
	double g = 0.5;
        uint idx = 0;
        bool save_to_file = false;
        double omega = 0.0;
        uint seed = 0;
        for (uint i = 1; i < argc; ++i) {
                if (strcmp(argv[i], "-g") == 0)
                        sscanf(argv[++i], "%lf", &g);
                else if (strcmp(argv[i], "-idx") == 0)
                        sscanf(argv[++i], "%u", &idx);
                else if (strcmp(argv[i], "-w") == 0)
                        sscanf(argv[++i], "%lf", &omega);
                else if (strcmp(argv[i], "-s") == 0)
                        save_to_file = true;
                else if (strcmp(argv[i], "-seed") == 0)
                        sscanf(argv[++i], "%u", &seed);
                else {
                        fprintf(stderr, "# ERROR: unknown argument %s.\n", argv[i]);
                        exit(EXIT_FAILURE);
                }
        }

        uint n_part = 6;
        uint tot_momentum = 6;
        FockStateSet set = fock_stateset_alloc_same_angular_momentum(tot_momentum, 130.0, 400, n_part);
        //FockStateSet set = fock_stateset_alloc(STATE_BESSEL, 200, n_part, 75.0);

        double *V, *T;
        if (save_to_file) {
                T = rgc_malloc(sizeof(*T)*set.size*set.size);
                V = rgc_malloc(sizeof(*V)*set.size*set.size);
                fock_matrix_fill(set, T, V);
                util_bin_tofile(T, sizeof(*T)*set.size*set.size, "T.dat");
                util_bin_tofile(V, sizeof(*V)*set.size*set.size, "V.dat");
                return 0;
        } else {
                V = util_bin_alloc_fromfile(sizeof(*V)*set.size*set.size, "V.dat");
                T = util_bin_alloc_fromfile(sizeof(*T)*set.size*set.size, "T.dat");
        }

        EigenSolutionSet sol = fock_eigen_vt_solve(g, set.size, T, V);
        FockCoeffStateSet cset = fock_coeffstateset_alloc_dcoeff(set, sol.vectors, 0.001);
        const uint nrandpos = 128;
        double* pos = fock_coeffstateset_alloc_random(cset, nrandpos);
        char fname[256];
        //sprintf(fname, "test.txt");
        sprintf(fname, "./results/data/random_diagonal/m6/rhoavg_g%.1lf.txt", g);
        FILE* f = rgc_fopen(fname, "w");
        for (double r = 0.0; r <= 1.0; r += 0.005) {
                fprintf(f, "%lf", r);
                double val[nrandpos];
                memset(val, 0, sizeof(*val)*nrandpos);
                for (double phi = 0.0; phi <= 2*M_PI; phi += 2*M_PI/60) {
                        double x = r*cos(phi);
                        double y = r*sin(phi);
                        double* psum = pos;
                        for (uint i = 0; i < nrandpos; ++i) {
                                psum[0] = x;
                                psum[6] = y;
                                double v = cabs(fock_coeffstateset_compute(cset, psum, psum+6));
                                val[i] += v*v;
                                psum += 2*n_part;
                        }
                }
                for (uint i = 0; i < nrandpos; ++i)
                        fprintf(f, " %e", val[i]);
                fputc('\n', f);
        }
        rgc_fclose(f);

        fock_stateset_free(&set);
	return 0;
}

