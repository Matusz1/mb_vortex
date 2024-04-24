#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include "src/fock.h"
#include <gsl/gsl_sf_bessel.h>
#include "src/gpe_radial.h"

int main(int argc, char **argv) {
        double omega = 0.0;
        //double mu, Epp;
        //BesselLinear blin = gpe_solve_m1(5.0*2.0, &mu, &Epp);
        //printf("radial: mu = %lf, E=%lf   vs   grid: %lf\n", 6.0*mu, 6.0*Epp, 58.109);
        //return 0;
        if (argc > 1)
                sscanf(argv[1], "%lf", &omega);
        //FILE* f_out = fopen("./results/data/energies/gpe-selfconsistent.txt", "w");
        uint size = 20;
        UVSol sol = bdg_solve(size, -2, 5*3.0, omega);
        for (uint i = 0; i < sol.size; ++i) {
                //printf("%15e    %15e\n", creal(sol.evals[i]), cimag(sol.evals[i]));
        }
        putchar('\n');
        printf("# ==================== PLUS FAMILY ======================== #\n");
        for (uint i = 0; i < 2*size; ++i) {
                if (!sol.plus_family[i])
                        continue;
                if (cimag(sol.evals[i]) == 0.0)
                        printf("%15lf\n", creal(sol.evals[i]));
                else
                        printf("(%15lf + %lfi)\n", creal(sol.evals[i]), cimag(sol.evals[i]));
                for (uint j = 0; j < sol.size; ++j) {
                        double v = cabs(sol.evecs[j*sol.size + i]);
                        v *= v;
                        BesselWaveState bws = fock_bessel_state_index(sol.types[j].nr);
                        if (v > 0.01) {
                                printf("\t%lf:   m:%d r:%d\n", v, bws.m, bws.zero_num);
                        }
                }
                putchar('\n');
        }
        putchar('\n');
        printf("# ==================== MINUS FAMILY ======================= #\n");
        for (uint i = 0; i < 2*size; ++i) {
                if (sol.plus_family[i])
                        continue;
                if (cimag(sol.evals[i]) == 0.0)
                        printf("%15lf\n", creal(sol.evals[i]));
                else
                        printf("(%15lf + %lfi)\n", creal(sol.evals[i]), cimag(sol.evals[i]));
                for (uint j = 0; j < sol.size; ++j) {
                        double v = cabs(sol.evecs[j*sol.size + i]);
                        v *= v;
                        BesselWaveState bws = fock_bessel_state_index(sol.types[j].nr);
                        if (v > 0.01) {
                                printf("\t%lf:   m:%d r:%d\n", v, bws.m, bws.zero_num);
                        }
                }
                putchar('\n');
        }
        //fclose(f_out);

	return 0;
}

