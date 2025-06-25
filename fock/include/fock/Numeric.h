#ifndef _NUMERIC_H_
#define _NUMERIC_H_

#include "Core.h"
#include <vector>

#define IS_EVEN(n) ((n)%2 == 0)
#define IS_ODD(n) ((n)%2 == 1)

namespace Mbs {

[[nodiscard]] std::vector<double> linspace(double start, double end, uint n, bool endpoint = true);
[[nodiscard]] std::vector<double> logspace(double start, double end, uint n);

[[nodiscard]] double simpson(double* f, uint n, double h);

inline void _dec_to_binarr(int* res, uint n, uint dim) {
    int pos = dim - 1;
    for (uint i = 0; i < dim+1; ++i)
        res[i] = 0;

    while (n != 0) {
        res[pos] = n & 1;
        res[dim] += res[pos];
        n /= 2;
        --pos;
    }
}

template<typename Ty>
[[nodiscard]] Ty permanent(Ty* A, uint n)
{
    int* chi = new int[n+1];
    const uint C = (1 << n); // 2^n
    Ty sum = 0;
    Ty rowsumprod, rowsum;

    for (uint k = 1; k < C; ++k) {
        rowsumprod = 1;
        _dec_to_binarr(chi, k, n);

        for (uint m = 0; m < n; ++m) {
            rowsum = 0;
            for (uint p = 0; p < n; ++p)
                rowsum += Ty(chi[p]) * A[m * n + p];
            rowsumprod *= rowsum;    
        }        

        sum += Ty((n-chi[n]) % 2 == 0 ? 1 : -1) * rowsumprod;
    }    

    delete[] chi;
    return sum;
}

} // namespace Mbs

#endif // _NUMERIC_H_
