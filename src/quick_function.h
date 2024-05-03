#ifndef _QUICK_FUNCTION_H_
#define _QUICK_FUNCTION_H_

/* QuickFunction - concept for quickly obtaining values
 * of bessel functions, interpolation and integraction
 * fucntions in range [0,1] */

#define QF_ARR_LEN 512
#define QF_DR (1.0/(QF_ARR_LEN-1))

typedef struct {
        double arr[QF_ARR_LEN];
} QuickFunction;

double qf_get_value(const QuickFunction* qf, double r);
double qf_integrate(const QuickFunction* qf);

#endif // !_QUICK_FUNCTION_H_
