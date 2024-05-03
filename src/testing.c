#include "testing.h"
#include "util.h"
#include "rho.h"

void testing_random(FockState fs) {
        printf("#| TEST: STARTING: testing_random\n"); fflush(stdout);
        uint npart = fs.basis.states[0].size;
        double p[2*npart];

        FILE* f_pos = util_fopen("test_random_pos.txt", "w");
        FILE* f_rho = util_fopen("test_random_rho.dat", "wb");

        /* Draw first positions */
        FockState fss[npart];
        Rho rhos[npart];
        fss[0] = fs;
        rhos[0] = rho_alloc(fs);

        rho_rand(rhos[0], FOCK_RHO_RAND_RELOAD, p, p + npart);


        printf("T| %u: fss.size = %u, rhos.basis.size = %u\n", 1, fss[0].basis.size, rhos[0].basis.size);

        for (uint j = 0; j < npart-1; ++j) {
                fss[j+1] = fock_operator_psi(fss[j], p[j], p[npart+j]);
                fock_coeffstateset_print(fss[j+1], 20);
                rhos[j+1] = rho_alloc(fss[j+1]);
                printf("T| %2u: fss.size = %u, rhos.basis.size = %u\n", j+2, fss[j+1].basis.size, rhos[j+1].basis.size);
                rho_rand(rhos[j+1], FOCK_RHO_RAND_RELOAD, &p[j+1], &p[npart+j+1]);
        }

        const uint nxy = 100;
        for (uint i = 0; i < npart; ++i) {
                fprintf(f_pos, "%lf %lf\n", p[i], p[npart+i]);
                for (int iy = 0; iy < nxy; ++iy) {
                        for (int ix = 0; ix < nxy; ++ix) {
                                double x = -1.0 + ix*(2.0)/(nxy-1);
                                double y = 1.0 - iy*(2.0)/(nxy-1);
                                double v = rho_compute_diag_xy(rhos[i], x, y);
                                util_fwrite(&v, 1, sizeof(v), f_rho);
                        }
                }
        }

        util_fclose(f_pos);
        util_fclose(f_rho);
        for (uint j = 1; j < npart; ++j)
                fock_state_free(fss+j);
        for (uint j = 0; j < npart; ++j)
                rho_free(rhos+j);

}

void testing_rho(FockState fs) {
        printf("# TEST: STARTING: testing_rho\n");
        Rho rho = rho_alloc(fs);

        FILE *f_out = util_fopen("test_rho.dat", "wb");

        uint nxy = 100;
        for (int iy = 0; iy < nxy; ++iy) {
                for (int ix = 0; ix < nxy; ++ix) {
                        double x = -1.0 + ix*(2.0)/(nxy-1);
                        double y = 1.0 - iy*(2.0)/(nxy-1);
                        double v = creal(rho_compute_xy(rho, x, y, x, y));
                        util_fwrite(&v, 1, sizeof(v), f_out);
                }
        }

        util_fclose(f_out);
        rho_free(&rho);
}
