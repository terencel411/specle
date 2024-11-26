#include "../globals.hpp"

void renorm(fftw_complex *data) {
    for (ptrdiff_t i = 0; i < local_n0; ++i) {
        for (ptrdiff_t j = 0; j < N1; ++j) {
            double *c = data[i*N1 + j];
            c[0] /= N0*N1;
            c[1] /= N0*N1;
        }
    }
}