#include "examples.h"
#include "rho.h"
#include "util.h"
#include "gpe_radial.h"

static double __search_cutoff(int m, uint npart, uint size_bound) {
        double cut = 10.0;
        while (true) {
                FockBasis basis = fock_basis_alloc_fixed_angular_momentum(npart, cut, npart);
                uint size = basis.size;
                fock_basis_free(&basis);
                if (size_bound < size)
                        break;
                cut += 1.0;
        }

        return cut - 1.0;
}

void example_random_sampling(double g, uint npart) {
        printf("# E| ===== STARTING EXAMPLE: %s\n", __func__);
        double cut = __search_cutoff(npart, npart, 3000);
        FockBasis basis = fock_basis_alloc_fixed_angular_momentum(npart, cut, npart);
        printf("# E| basis.size = %u\n", basis.size);
        printf("# E| Diagonalazing ...\n");
        DEigenSol esol = fock_eigen_solve(g, basis);
        FockState state = {
                .basis = basis,
                .a = esol.U,
        };

        printf("# E| Drawing positions ...\n");
        double* pos = fock_state_alloc_random(state, 8);

        /* Save conditional wave functions to output files */
        printf("# E| Saving wave functions ...\n");
        FILE* f_wf = util_fopen("example_output.dat", "wb");
        const uint nxy = 100;
        complex double out[nxy*nxy];
        for (uint i = 0; i < 8; ++i) {
                printf(" [ %2u / %2u ]\n", i+1, 8);
                uint ixy = 0;
                for (int iy = 0; iy < nxy; ++iy) {
                        for (int ix = 0; ix < nxy; ++ix) {
                                double x = -1.0 + ix*(2.0)/(nxy-1);
                                double y = 1.0 - iy*(2.0)/(nxy-1);
                                pos[i*(2*npart)] = x;
                                pos[i*(2*npart) + npart] = y;
                                out[ixy] = fock_state_compute(state, &pos[i*(2*npart)], &pos[i*(2*npart)+npart]);
                                ++ixy;
                        }
                }
                util_fwrite(out, nxy*nxy, sizeof(*out), f_wf);
        }
        util_fclose(f_wf);

        printf("# E| Finishing ...\n");
        free(pos);
        fock_basis_free(&basis);
}

void example_gpe() {

}

void example_bdg(double omega) {
        //double mu, Epp;
        //BesselLinear blin = gpe_solve_m1(5.0*2.0, &mu, &Epp);
        //printf("radial: mu = %lf, E=%lf   vs   grid: %lf\n", 6.0*mu, 6.0*Epp, 58.109);
        //return 0;
        //FILE* f_out = fopen("./results/data/energies/gpe-selfconsistent.txt", "w");
        uint size = 20;
        UVSol sol = bdg_solve(size, 0, 5*3.0, omega);
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
                        StateInfo bws = state_info(sol.types[j].nr);
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
                        StateInfo bws = state_info(sol.types[j].nr);
                        if (v > 0.01) {
                                printf("\t%lf:   m:%d r:%d\n", v, bws.m, bws.zero_num);
                        }
                }
                putchar('\n');
        }
}
