#include <tgmath.h>

#include "examples.h"
#include "rho.h"
#include "util.h"
#include "gpe_radial.h"

static double __search_cutoff(int m, uint npart, uint size_bound) {
        double cut = 10.0;
        while (true) {
                FockBasis basis = fock_basis_alloc_fixed_angular_momentum(m, cut, npart);
                uint size = basis.size;
                fock_basis_free(&basis);
                if (size_bound < size)
                        break;
                cut += 1.0;
        }

        return cut - 1.0;
}

void example_random_sampling(double gamma, uint npart) {
        printf("# E| ===== STARTING EXAMPLE: %s\n", __func__);
        int m = 6*npart;
        double cut = __search_cutoff(m, npart, 3500);
        FockBasis basis = fock_basis_alloc_fixed_angular_momentum(m, cut, npart);
        printf("# E| basis.size = %u\n", basis.size);
        printf("# E| Diagonalazing ...\n");
        double g = gamma / npart;
        DEigenSol esol = fock_eigen_solve(g, basis);
        FockState state = {
                .basis = basis,
                .a = esol.U,
        };

        printf("# E| Drawing positions ...\n");
        double* pos = fock_state_alloc_random(state, 8);

        /* Save conditional wave functions to output files */
        printf("# E| Saving wave functions ...\n");
        char fname[256];
        sprintf(fname, "draws_n%u_gN%.1lf_m%d.dat", npart, gamma, m);
        FILE* f_wf = util_fopen(fname, "wb");
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

void example_gpe_energy(uint npart) {
        char fname[256];
        sprintf(fname, "gpe_Epp_n%u.txt", npart);
        FILE* f_out = util_fopen(fname, "w");
        fprintf(f_out, "# gN       mu        Epp\n");
        for (double gamma = 0.0; gamma < 18.01; gamma += 0.1) {
                double g = gamma/npart;
                double mu, Epp;
                gpe_solve_m0(g*(npart-1), &mu, &Epp);
                fprintf(f_out, "%lf %lf %lf\n", gamma, mu, Epp);
        }
        util_fclose(f_out);
}

void example_bdg_write(int m, double gN, bool vortex_state) {
        uint size = 20;
        UVSol sol;
        if (vortex_state) 
                sol = bdg_solve(size, m, gN);
        else
                sol = bdg_solve_m0(size, m, gN);
        printf("# ====================== PLUS FAMILY ====================== #\n");
        for (uint i = 0; i < 2*size; ++i) {
                if (!sol.plus_family[i])
                        continue;

                if (cimag(sol.evals[i]) == 0.0)
                        printf("%15lf\n", creal(sol.evals[i]));
                else
                        printf("(%15lf + %lfi)\n", creal(sol.evals[i]), cimag(sol.evals[i]));

                for (uint j = 0; j < sol.size; ++j) {
                        double v = cabs(sol.evecs[i*sol.size + j]);
                        v *= v;
                        if (v > 0.001) {
                                StateInfo si = state_info(sol.types[j].nr);
                                printf("\t%lf:   m:%d r:%d\n", v, si.m, si.zero_num);
                        }
                }
                putchar('\n');
        }
        putchar('\n');
        printf("# ====================== MINUS FAMILY ====================== #\n");
        for (uint i = 0; i < 2*size; ++i) {
                if (sol.plus_family[i])
                        continue;

                if (cimag(sol.evals[i]) == 0.0)
                        printf("%15lf\n", creal(sol.evals[i]));
                else
                        printf("(%15lf + %lfi)\n", creal(sol.evals[i]), cimag(sol.evals[i]));

                for (uint j = 0; j < sol.size; ++j) {
                        double v = cabs(sol.evecs[i*sol.size + j]);
                        v *= v;
                        if (v > 0.01) {
                                StateInfo si = state_info(sol.types[j].nr);
                                printf("\t%lf:   m:%d r:%d\n", v, si.m, si.zero_num);
                        }
                }
                putchar('\n');
        }
}

void example_bdg(uint npart, double gamma_max, uint ncuts) {
        char fname[256];
        sprintf(fname, "bdg_n%u.txt", npart);
        FILE *f_bdg = util_fopen(fname, "w");
        fprintf(f_bdg, "# gamma  lambda0  E_GPE    E_BdG\n");

        for (uint ig = 0; ig < ncuts; ++ig) {
                // Determine interaction strength
                double gamma = ig * gamma_max/(ncuts-1);
                double g = gamma / npart;
                uint size = 20;
                double depletion = 0.0;
                double e_corr = 0.0;
                double mu, Epp;
                gpe_solve_m1(g*(npart-1), &mu, &Epp);
                for (int m = -6; m <= 4; ++m) {
                        UVSol sol = bdg_solve(size, m, (npart-1)*g);
                        for (uint i = 0; i < 2*size; ++i) {
                                if (!sol.plus_family[i])
                                        continue;
                                if (cimag(sol.evals[i]) != 0.0)
                                        continue;
                                for (uint j = 0; j < sol.size; ++j) {
                                        if (!sol.types[j].isUp) {
                                                double v = cabs(sol.evecs[i*sol.size + j]);
                                                e_corr -= creal(sol.evals[i])*v*v;
                                                if (cimag(sol.evals[i]) == 0.0)
                                                        depletion += v*v;
                                        }
                                }
                        }
                        uvsol_free(&sol);
                }
                fprintf(f_bdg, "%lf %lf %lf %lf\n", gamma, 1.0 - (depletion/npart), npart*Epp, npart*Epp+e_corr);
        }
        util_fclose(f_bdg);
}

void example_rho_diag(uint npart, double gamma_min, double gamma_max, uint ncuts, uint ssize) {
        int m = 0;
        double cut = __search_cutoff(m, npart, ssize);
        FockBasis basis = fock_basis_alloc_fixed_angular_momentum(m, cut, npart);

        printf("# size = %u\n", basis.size);
        double *T = util_malloc(sizeof(*T)*basis.size*basis.size);
        double *V = util_malloc(sizeof(*V)*basis.size*basis.size);
        fock_matrix_fill(basis, T, V);

        char fname[256];
        sprintf(fname, "mb_depletion_n%u_m%u.txt", npart, m);
        FILE *f_rho = util_fopen(fname, "a");

        for (uint i = 0; i < ncuts; ++i) {
                // Determine interaction strength
                double gamma = gamma_min + i * (gamma_max-gamma_min)/(ncuts-1);
                double g = gamma / npart;

                // Solve eigenproblem
                DEigenSol fsol = fock_eigen_vt_solve(g, basis.size, T, V);

                // Diagonalize density matrix
                Rho rho = rho_alloc((FockState){.basis = basis, .a = fsol.U});
                DEigenSol rsol = rho_diag(&rho);

                // Determine depletion
                double lambda0 = rsol.val[rsol.size - 1]/npart;
                double epp = fsol.val[0];
                printf("gamma = %lf, g = %lf, lambda0 = %lf\n", gamma, g, lambda0);
                fprintf(f_rho, "%lf %lf %lf\n", gamma, epp/npart, lambda0);
                fflush(f_rho);

                rho_free(&rho);
                deigen_sol_free(&fsol);
                deigen_sol_free(&rsol);
        }

        util_fclose(f_rho);
        util_free(T);
        util_free(V);
        fock_basis_free(&basis);
}

void example_rho_evec(uint npart, double gamma, uint ssize) {
        int m = npart;
        double cut = __search_cutoff(m, npart, ssize);
        FockBasis basis = fock_basis_alloc_fixed_angular_momentum(m, cut, npart);

        DEigenSol fsol = fock_eigen_solve(gamma/npart, basis);
        Rho rho = rho_alloc((FockState){.basis = basis, .a = fsol.U});
        DEigenSol rsol = rho_diag(&rho);

        char fname[256];
        sprintf(fname, "rho_vec_%u_%.1lf.txt", npart, gamma);
        FILE* f_out = util_fopen(fname, "w");

        for (uint i = 0; i < 6; ++i) {
                fprintf(f_out, "%lf\t", rsol.val[rsol.size - 1 - i]/npart);
                for (double r = -1.0; r < 1.001; r += 0.01) {
                        complex double v = 0.0;
                        for (uint j = 0; j < rsol.size; ++j) {
                                complex double a = rsol.U[(rsol.size-1-i)*rsol.size + j];
                                v += a*state_value_r(rho.basis.states[j], fabs(r));
                        }
                        fprintf(f_out, " %e", cabs(v*v));
                }
                fputc('\n', f_out);
        }

        util_fclose(f_out);

        rho_free(&rho);
        deigen_sol_free(&fsol);
        deigen_sol_free(&rsol);
}

void example_rho_gpe_compariton(uint npart, double gamma, uint ssize) {
        int m = npart/2;
        double cut = __search_cutoff(m, npart, ssize);
        FockBasis basis = fock_basis_alloc_fixed_angular_momentum(m, cut, npart);

        printf("# E| Diagonalizing ...\n");
        double g = gamma/npart;
        DEigenSol fsol = fock_eigen_solve(g, basis);
        Rho rho = rho_alloc((FockState){.basis = basis, .a = fsol.U});
        fock_basis_free(&basis);
        deigen_sol_free(&fsol);

        printf("# E| Solving GPE ...\n");
        double mu, Epp;
        BesselLinear bl = gpe_solve_m1(g*(npart-1), &mu, &Epp);
        QuickFunction gpe_qf = gpe_bes_lin_get_qf(&bl);

        printf("# E| Saving to file ...\n");
        char fname[256];
        sprintf(fname, "mb_gpe_comparison_n%u_m%d_gN%.1lf.txt", npart, m, gamma);
        FILE *f_out = util_fopen(fname, "w");
        fprintf(f_out, "# MB     GPE\n");
        for (double r = -1.0; r <= 1.001; r += 0.01) {
                double v_mb = rho_compute_diag(rho, fabs(r), 0.0)/npart;
                double v_gpe = qf_get_value(&gpe_qf, fabs(r));
                v_gpe *= v_gpe;
                fprintf(f_out, "%lf %lf\n", v_mb, v_gpe);
        }
        util_fclose(f_out);

        rho_free(&rho);
}

void example_bdg_uv(uint npart, double gamma) {
        double g = gamma/npart;
        uint size = 20;
        double mu, Epp;
        BesselLinear bl = gpe_solve_m1(g*(npart-1), &mu, &Epp);
        QuickFunction bqf = gpe_bes_lin_get_qf(&bl);
        UVSol sol = bdg_solve(size, -2, g*(npart-1));

        char fname[256];
        sprintf(fname, "uv_state.txt");
        FILE *f_out = util_fopen(fname, "w");
        for (uint i = 0; i < 2*size; ++i) {
                if (!sol.plus_family[i])
                        continue;
                if (cimag(sol.evals[i]) != 0.0)
                        continue;
                if (creal(sol.evals[i]) > 0.0)
                        continue;

                printf("%lf\n", creal(sol.evals[i]));
                for (uint j = 0; j < sol.size; ++j) {
                        double v = cabs(sol.evecs[i*sol.size + j]);
                        v *= v;
                        if (v > 0.001) {
                                StateInfo si = state_info(sol.types[j].nr);
                                printf("\t%lf:   m:%d r:%d\n", v, si.m, si.zero_num);
                        }
                }

                for (double r = -1.0; r <= 1.001; r += 0.01) {
                        complex double tot = 0.0;
                        for (uint j = 0; j < sol.size; ++j) {
                                double a = creal(sol.evecs[i*sol.size + j]);
                                double v = a*state_value_r(sol.types[j].nr, fabs(r));
                                tot += v;
                        }
                        double v = fabs(tot);
                        v *= v;
                        double gpev = qf_get_value(&bqf, fabs(r));
                        gpev *= gpev;
                        fprintf(f_out, "%lf %lf %lf\n", r, v, gpev);
                }
        }
        util_fclose(f_out);
        uvsol_free(&sol);
}
