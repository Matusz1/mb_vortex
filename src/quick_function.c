#include "quick_function.h"
#include "elementary.h"
#include <tgmath.h>

double qf_get_value(const QuickFunction* qf, double r) {
        r = fabs(r);
        if (r > 1.0)
                return 0.0;
        uint l = (uint) (r / QF_DR);
        uint u = l + 1;
        double vl = qf->arr[l];
        double vu = qf->arr[u];
        double dfdr = (vu - vl) / QF_DR;
        double diff = (r - l*QF_DR) * dfdr;
        return vl + diff;
}

double qf_integrate(const QuickFunction* qf) {
        double sum = 0.0;
        for (uint i = 1; i < QF_ARR_LEN-1; ++i) {
                double r = i*QF_DR;
                sum += r*qf->arr[i];
        }
        sum += 0.5*qf->arr[QF_ARR_LEN-1];
        return sum*QF_DR;
}
