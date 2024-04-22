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

        G1State gfunc = fock_g1_alloc(set, sol.vectors);

        // Print G function
        uint ixy = 0;
        double rho_arr[256*256];
        for (uint iy = 0; iy < 256; ++iy) {
                for (uint ix = 0; ix < 256; ++ix) {
                        double y = -1.28 + iy*(0.01);
                        double x = -1.28 + ix*(0.01);
                        y = ((double)iy - 128.0) / 100.0;
                        x = ((double)ix - 128.0) / 100.0;
                        rho_arr[ixy] = cabs(fock_g1_compute(&gfunc, x, y, x, y));
                        ++ixy;
                }
        }
        char fname_rho[256];
        sprintf(fname_rho, "./results/data/MB_rho/rho_n%u_m%u_g%.1lf.dat", n_part, tot_momentum, g);
        util_bin_tofile(rho_arr, 256*256*sizeof(*rho_arr), fname_rho);
        return 0;

        // Diagonalize G function and print top 5 results
        uint nx = 40;
        uint ny = 40;
        Lattice2D lattice = fock_lattice_alloc((Vec2d){.x=-1.0,.y=-1.0}, (Vec2d){.x=1.0,.y=1.0}, nx, ny);
        complex double* mat = fock_g1_matrix_alloc(lattice, &gfunc);
        CmplxEigenSolutionSet gsol = fock_g1_matrix_solve(mat, lattice.size);

        char fname_diag[256];
        sprintf(fname_diag, "./results/data/p6m%u/diag_g%.1lf.dat", tot_momentum, g);
        FILE* f_diag = rgc_fopen(fname_diag, "wb");
        rgc_fwrite(gsol.values+(gsol.size-5), sizeof(*gsol.values), 5, f_diag);
        rgc_fwrite(gsol.vectors+gsol.size*(gsol.size-5), sizeof(*gsol.vectors), 5*gsol.size, f_diag);
        rgc_fclose(f_diag);
        

        fock_stateset_free(&set);
	return 0;
}

