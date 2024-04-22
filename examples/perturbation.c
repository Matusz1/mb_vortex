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
        for (uint i = 1; i < argc; ++i) {
                if (strcmp(argv[i], "-g") == 0)
                        sscanf(argv[++i], "%lf", &g);
                else if (strcmp(argv[i], "-idx") == 0)
                        sscanf(argv[++i], "%u", &idx);
                else if (strcmp(argv[i], "-w") == 0)
                        sscanf(argv[++i], "%lf", &omega);
                else if (strcmp(argv[i], "-seed") == 0)
                        sscanf(argv[++i], "%u", &seed);
                else {
                        fprintf(stderr, "# ERROR: unknown argument %s.\n", argv[i]);
                        exit(EXIT_FAILURE);
                }
        }

        uint n_part = 6;
        FockStateSet set0 = fock_stateset_alloc_same_angular_momentum(0, 160.0, 400, n_part);
        FockStateSet set1 = fock_stateset_alloc_same_angular_momentum(n_part, 170.0, 400, n_part);
        //printf("%u %u\n", set0.size, set1.size);
        //return 0;

        FockState state0 = fock_init_quick_bessel(6, (int[6]){0,0,0,0,0,0}, (uint[6]){1,1,1,1,1,1});
        FockState state1 = fock_init_quick_bessel(6, (int[6]){1,1,1,1,1,1}, (uint[6]){1,1,1,1,1,1});

        double e0 = fock_energy(&state0);
        double e1 = fock_energy(&state1);
        fprintf(stderr, "OK1\n");
        double v01 = fock_stateset_perturbation_energy_one(1.0, &state0);
        double v02 = fock_stateset_perturbation_energy_two(1.0, set0, &state0);
        double v03 = fock_stateset_perturbation_energy_three(1.0, set0, &state0);
        fprintf(stderr, "OK2\n");
        double v11 = fock_stateset_perturbation_energy_one(1.0, &state1);
        double v12 = fock_stateset_perturbation_energy_two(1.0, set1, &state1);
        double v13 = fock_stateset_perturbation_energy_three(1.0, set1, &state1);
        fprintf(stderr, "OK3\n");
        FILE* f = rgc_fopen("./results/data/perturbation.txt", "w");
        for (double gg = 0.0; gg <= 3.00; gg += 0.01) {
                fprintf(f, "%.2lf %lf %lf %lf %lf %lf %lf\n", gg, e0+v01*gg, e0+v01*gg+v02*gg*gg, e0+v01*gg+v02*gg*gg+v03*gg*gg*gg,
                        e1+v11*gg, e1+v11*gg+v12*gg*gg, e1+v11*gg+v12*gg*gg+v13*gg*gg*gg);
                fflush(f);
        }
        rgc_fclose(f);

        fock_stateset_free(&set0);
        fock_stateset_free(&set1);
	return 0;
}

