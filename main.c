#include "src/util.h"
#include "src/fock.h"
#include <math.h>

int main(int argc, char **argv) {

        for (uint i = 0; i < 200; ++i) {
                double ener = fock_single_energy(i, STATE_BESSEL);
                printf("%4u: %10lf\n", i, ener);
        }
        return 0;
        uint npart = 6;
        FockStateSet set = fock_stateset_alloc_same_angular_momentum(6, 100.0, 400, npart);
        printf("# set.size = %u\n", set.size);
        // return 0;
        EigenSolutionSet sol = fock_eigen_solve(3.0, set);
        FockCoeffStateSet cset = fock_coeffstateset_alloc_dcoeff(set, sol.vectors, 1.0e-4);
        uint ndraws = 16;
        double* pos = fock_coeffstateset_alloc_random(cset, ndraws);

        for (uint i = 0; i < ndraws; ++i) {
                for (uint j = 0; j < npart; ++j) {
                        double x = pos[i*(2*npart)+j];
                        double y = pos[i*(2*npart)+j+npart];
                        printf("%lf %lf\n", x, y);
                }
        }

        free(pos);
        fock_coeffstateset_free(&cset);
        fock_stateset_free(&set);
	return 0;
}

