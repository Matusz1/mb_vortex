#include "src/util.h"
#include "src/fock.h"
#include "src/rho.h"
#include "src/testing.h"
#include "src/examples.h"
#include <math.h>

int main(int argc, char **argv) {
        //example_random_sampling(2.0, 6);
        double omega = 0.0;
        if (argc > 1)
                sscanf(argv[1], "%lf", &omega);
        example_bdg(omega);

	return 0;
}

