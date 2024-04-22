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
        double omega = 0.0;
        uint seed = 0;
        bool check_size = false;
        uint tot_momentum = 0;
        double Ecut = 100.0;
        for (uint i = 1; i < argc; ++i) {
                if (strcmp(argv[i], "-g") == 0)
                        sscanf(argv[++i], "%lf", &g);
                else if (strcmp(argv[i], "-idx") == 0)
                        sscanf(argv[++i], "%u", &idx);
                else if (strcmp(argv[i], "-w") == 0)
                        sscanf(argv[++i], "%lf", &omega);
                else if (strcmp(argv[i], "-c") == 0)
                        check_size = true;
                else if (strcmp(argv[i], "-seed") == 0)
                        sscanf(argv[++i], "%u", &seed);
                else if (strcmp(argv[i], "-m") == 0)
                        sscanf(argv[++i], "%u", &tot_momentum);
                else if (strcmp(argv[i], "-e") == 0)
                        sscanf(argv[++i], "%lf", &Ecut);
                else {
                        fprintf(stderr, "# ERROR: unknown argument %s.\n", argv[i]);
                        exit(EXIT_FAILURE);
                }
        }

        uint n_part = 6;
        FockStateSet set = fock_stateset_alloc_same_angular_momentum(tot_momentum, Ecut, 400, n_part);
        if (check_size) {
                printf("%u\n", set.size);
                return 0;
        }
        //FockStateSet set = fock_stateset_alloc(STATE_BESSEL, 200, n_part, 75.0);

        uint saved_npart;
        uint saved_mom;
        uint saved_size;
        FILE* file_saved_data = rgc_fopen("./TV.txt", "r");
        rgc_error(fscanf(file_saved_data, "%u %u %u", &saved_npart, &saved_mom, &saved_size) != 3, "TV.txt\n");
        rgc_fclose(file_saved_data);
        double *V, *T;
        if (saved_size != set.size || saved_npart != n_part || saved_mom != tot_momentum) {
                printf("# Saving T and V matrices to file ...\n");
                T = rgc_malloc(sizeof(*T)*set.size*set.size);
                V = rgc_malloc(sizeof(*V)*set.size*set.size);
                fock_matrix_fill(set, T, V);
                util_bin_tofile(T, sizeof(*T)*set.size*set.size, "T.dat");
                util_bin_tofile(V, sizeof(*V)*set.size*set.size, "V.dat");

                FILE* file_saved_data = rgc_fopen("./TV.txt", "w");
                fprintf(file_saved_data, "%u %u %u\n", n_part, tot_momentum, set.size);
                rgc_fclose(file_saved_data);

                printf("# Saving matrices done.\n");
        } else {
                printf("# Loading T and V matrices from file ...\n");
                V = util_bin_alloc_fromfile(sizeof(*V)*set.size*set.size, "V.dat");
                T = util_bin_alloc_fromfile(sizeof(*T)*set.size*set.size, "T.dat");
                printf("# Loading matrices done.\n");
        }

        EigenSolutionSet sol = fock_eigen_vt_solve(g, set.size, T, V);
        FockCoeffStateSet cset = fock_coeffstateset_alloc_dcoeff(set, sol.vectors, 1.0e-4);
        double* pos = fock_coeffstateset_alloc_random(cset, 8);

        complex double output[8*100*100];
        uint ixy = 0;
        for (uint i = 0; i < 8; ++i) {
                for (uint iy = 0; iy < 100; ++iy) {
                        for (uint ix = 0; ix < 100; ++ix) {
                                double x = -1.0 + 2.0*ix/99.0;
                                double y =  1.0 - 2.0*iy/99.0;
                                pos[0] = x;
                                pos[n_part] = y;
                                output[ixy] = fock_coeffstateset_compute(cset, pos, pos+n_part);
                                ++ixy;
                        }
                }
                pos += n_part*2;
        }
        char fname_output[256];
        sprintf(fname_output, "./results/data/draws/draw_g%.1lf_m%u.dat", g, tot_momentum);
        util_bin_tofile(output, sizeof(output), fname_output);
        

        fock_stateset_free(&set);
	return 0;
}

