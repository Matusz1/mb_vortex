#ifndef _IMPORTANCE_TRUNCATION_H_
#define _IMPORTANCE_TRUNCATION_H_

#include "Core.h"
#include "Workspace.h"

namespace Mbs {

struct TruncationParameters {
    double kappa_min{1e-5};
    double C_min{std::numeric_limits<double>::min()};
    double ener_cutoff{100.0};
};

Workspace importance_truncation_scheme(
    double g,
    uint N,
    int mom,
    TruncationParameters params
);

}

#endif // _IMPORTANCE_TRUNCATION_H_
